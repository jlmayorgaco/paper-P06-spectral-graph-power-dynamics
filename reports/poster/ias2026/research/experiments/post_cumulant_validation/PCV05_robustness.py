# ruff: noqa: E501  -- table labels kept on one line
"""PCV05 - Phase 7 targeted robustness over deterministic stress-test envelopes.

Preregistration 5d0b1986, section 9. Coverage fractions over declared envelopes are
NOT probabilities.

- EM-f: E37 machine envelope, one factor per group applied to every SG.
- EM-u: same bounds, an independent factor per SG per group.
- EC: converter stress test (no project source), an independent factor per converter
  per group; Q/V gains excluded (policy coordinate).
- EMC: EM-u together with EC.

Latin hypercube N = 100 per envelope (seed 20260917 + envelope index) plus
one-at-a-time extremes. Policy fixed at P4. Per draw: the 16 subsets of H4 through the
FC01 direct path, H, kappa, alpha(H4), alpha(H4, g = 1) and g* by bisection on
[0.03625, 1.0] (25 steps) when H4 is unstable at P4 and stable at g = 1.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _pcv import (
    H4,
    POINTS,
    SUBSETS4,
    hypergraph,
    label,
    out_dir,
    transverse_of,
    write_json,
)
from scipy.stats import qmc

from _f7_common import LEAK
from _overnight import pin_blas_threads
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

OUT = out_dir("PCV05")
P4 = POINTS["P4"]
SG_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
CONV_BUSES = H4
N_LHS = 100
SEED = 20260917
G_LO, G_HI, G_STEPS = 0.03625, 1.0, 25
MACHINE_GROUPS = {
    "M": (("m",), 0.80, 1.20),
    "XP": (("xd1", "xq1"), 0.90, 1.10),
    "KA": (("ka",), 0.85, 1.15),
    "TE": (("ta",), 0.85, 1.15),
    "PSS_K": (("pss_gain",), 0.80, 1.20),
    "PSS_T": (("pss_washout", "pss_wash_lag", "pss_lag"), 0.80, 1.20),
}
CONV_GROUPS = {
    "PLL": (("kp_pll", "ki_pll"), 0.80, 1.20),
    "OUTER": (("kp_p", "ki_p", "kp_q", "ki_q"), 0.80, 1.20),
    "CURRENT": (("kp_i", "ki_i"), 0.80, 1.20),
    "TAU_P": (("tau_p",), 0.80, 1.20),
    "XF": (("xf",), 0.90, 1.10),
}
ENVELOPES = ("EM-f", "EM-u", "EC", "EMC")
DEFAULT_CONV = ConverterParameters()


# ------------------------------------------------------------------ draws --
def lhs(dim, index):
    return qmc.LatinHypercube(d=dim, rng=np.random.default_rng(SEED + index)).random(
        N_LHS
    )


def scale(u, lo, hi):
    return float(lo + u * (hi - lo))


def make_draws():
    draws = [{"envelope": "NOMINAL", "draw": 0, "fleet": {}, "unit": {}, "conv": {}}]
    mg, cg = list(MACHINE_GROUPS), list(CONV_GROUPS)
    # EM-f
    for i, row in enumerate(lhs(len(mg), 0)):
        draws.append(
            {
                "envelope": "EM-f",
                "draw": i,
                "unit": {},
                "conv": {},
                "fleet": {
                    g: scale(row[j], *MACHINE_GROUPS[g][1:]) for j, g in enumerate(mg)
                },
            }
        )
    # EM-u, EC, EMC
    nu, nc = len(SG_BUSES) * len(mg), len(CONV_BUSES) * len(cg)

    def unit_of(row):
        return {
            b: {
                g: scale(row[i * len(mg) + j], *MACHINE_GROUPS[g][1:])
                for j, g in enumerate(mg)
            }
            for i, b in enumerate(SG_BUSES)
        }

    def conv_of(row):
        return {
            b: {
                g: scale(row[i * len(cg) + j], *CONV_GROUPS[g][1:])
                for j, g in enumerate(cg)
            }
            for i, b in enumerate(CONV_BUSES)
        }

    for i, row in enumerate(lhs(nu, 1)):
        draws.append(
            {
                "envelope": "EM-u",
                "draw": i,
                "fleet": {},
                "conv": {},
                "unit": unit_of(row),
            }
        )
    for i, row in enumerate(lhs(nc, 2)):
        draws.append(
            {"envelope": "EC", "draw": i, "fleet": {}, "unit": {}, "conv": conv_of(row)}
        )
    for i, row in enumerate(lhs(nu + nc, 3)):
        draws.append(
            {
                "envelope": "EMC",
                "draw": i,
                "fleet": {},
                "unit": unit_of(row[:nu]),
                "conv": conv_of(row[nu:]),
            }
        )
    # one-at-a-time extremes
    k = 0
    for g, (_, lo, hi) in MACHINE_GROUPS.items():
        for v in (lo, hi):
            draws.append(
                {
                    "envelope": "OAT-M",
                    "draw": k,
                    "fleet": {g: v},
                    "unit": {},
                    "conv": {},
                    "oat": f"{g}={v}",
                }
            )
            k += 1
    k = 0
    for g, (_, lo, hi) in CONV_GROUPS.items():
        for v in (lo, hi):
            draws.append(
                {
                    "envelope": "OAT-C",
                    "draw": k,
                    "fleet": {},
                    "unit": {},
                    "conv": {b: {g: v} for b in CONV_BUSES},
                    "oat": f"{g}={v}",
                }
            )
            k += 1
    return draws


def case_kwargs(draw, theta):
    g, k, t, h = theta
    assert h == 1.0
    ms = {"ka": k, "ta": t}
    for grp, f in draw["fleet"].items():
        for field in MACHINE_GROUPS[grp][0]:
            ms[field] = ms.get(field, 1.0) * f
    mbs = None
    if draw["unit"]:
        mbs = {
            b: {
                field: f
                for grp, f in groups.items()
                for field in MACHINE_GROUPS[grp][0]
            }
            for b, groups in draw["unit"].items()
        }
    base = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK)
    convs = None
    if draw["conv"]:
        convs = {
            b: replace(
                base,
                **{
                    field: getattr(DEFAULT_CONV, field) * f
                    for grp, f in groups.items()
                    for field in CONV_GROUPS[grp][0]
                },
            )
            for b, groups in draw["conv"].items()
        }
    return {
        "converter": base,
        "converters": convs,
        "machine_scaling": ms,
        "machine_bus_scaling": mbs,
    }


def solve(draw, members, theta):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}), **case_kwargs(draw, theta)
    )


# ----------------------------------------------------------------- worker --
def evaluate(draw):
    pin_blas_threads()
    rec = {
        "envelope": draw["envelope"],
        "draw": draw["draw"],
        "oat": draw.get("oat", ""),
        "factors": json.dumps(
            {"fleet": draw["fleet"], "unit": draw["unit"], "conv": draw["conv"]},
            sort_keys=True,
        ),
    }
    statuses = {}
    for s in SUBSETS4:
        try:
            tr = transverse_of(solve(draw, s, P4))
        except Exception as exc:  # infeasible power flow / initialization
            tr = {
                "status": "SOLVE_FAILED",
                "alpha": float("nan"),
                "crit_hz": float("nan"),
                "error": repr(exc),
            }
        statuses[s] = tr["status"]
        rec[f"alpha_{label(s)}"] = tr["alpha"]
        rec[f"status_{label(s)}"] = tr["status"]
        if s == H4:
            rec["alpha_H4"], rec["crit_hz_H4"] = tr["alpha"], tr["crit_hz"]
    failed = any(v == "SOLVE_FAILED" for v in statuses.values())
    hg = hypergraph(
        {
            s: ("BOUNDARY_OR_UNRESOLVED" if v == "SOLVE_FAILED" else v)
            for s, v in statuses.items()
        }
    )
    rec.update(
        {"H": hg["H"], "kappa": hg["kappa"], "exact": hg["exact"] and not failed}
    )
    singles_stable = all(statuses[(b,)] == "STABLE" for b in H4)
    base_stable = statuses[()] == "STABLE"
    rec["base_stable"] = base_stable
    rec["non_composable"] = bool(
        base_stable and singles_stable and hg["kappa"] in (2.0, 3.0, 4.0)
    )
    rec["H_is_H4"] = hg["H"] == label(H4)
    rec["kappa_class"] = (
        "base-unstable"
        if not base_stable
        else ("inf" if not np.isfinite(hg["kappa"]) else str(int(hg["kappa"])))
    )

    def a_h4(g):
        return transverse_of(solve(draw, H4, (g, P4[1], P4[2], P4[3])), classify=False)[
            "alpha"
        ]

    a_hi = a_h4(G_HI)
    rec["alpha_H4_g1"] = a_hi
    a_lo = rec["alpha_H4"]
    if a_lo > 0 > a_hi:
        lo, hi = G_LO, G_HI
        for _ in range(G_STEPS):
            mid = 0.5 * (lo + hi)
            if a_h4(mid) > 0:
                lo = mid
            else:
                hi = mid
        rec["g_star"] = 0.5 * (lo + hi)
        rec["crossing"] = "crossing"
    else:
        rec["g_star"] = float("nan")
        rec["crossing"] = "no crossing"
    return rec


def decision(x):
    return "PASS" if x >= 0.8 else ("PARTIAL" if x >= 0.5 else "FAIL")


def main(argv) -> int:
    started = time.time()
    draws = make_draws()
    with Pool(16) as pool:
        rows = pool.map(evaluate, draws, chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PCV05_draws.csv", index=False)
    nominal = df[df.envelope == "NOMINAL"].iloc[0]
    summ = []
    for env in ("NOMINAL", *ENVELOPES, "OAT-M", "OAT-C"):
        d = df[df.envelope == env]
        n = len(d)
        r1 = float(d.non_composable.mean())
        r2 = float(d.H_is_H4.mean())
        kd = d.kappa_class.value_counts().to_dict()
        gs = d.g_star.dropna()
        direction = float((d.alpha_H4_g1 < d.alpha_H4).mean())
        summ.append(
            {
                "envelope": env,
                "n_draws": n,
                "n_unresolved_or_failed": int((~d.exact.astype(bool)).sum()),
                "R1_non_composable_fraction": r1,
                "R1_decision": decision(r1) if env in ENVELOPES else "",
                "R2_H_equals_H4_fraction": r2,
                "R2_decision": decision(r2) if env in ENVELOPES else "",
                "R3_kappa_distribution": json.dumps(
                    {k: int(v) for k, v in sorted(kd.items())}
                ),
                "R4_n_crossing": int(len(gs)),
                "R4_g_star_median": float(gs.median()) if len(gs) else float("nan"),
                "R4_g_star_q25": float(gs.quantile(0.25)) if len(gs) else float("nan"),
                "R4_g_star_q75": float(gs.quantile(0.75)) if len(gs) else float("nan"),
                "R4_g_star_min": float(gs.min()) if len(gs) else float("nan"),
                "R4_g_star_max": float(gs.max()) if len(gs) else float("nan"),
                "R4_direction_fraction": direction,
                "R4_direction_robust": bool(direction >= 0.9),
                "alpha_H4_min": float(d.alpha_H4.min()),
                "alpha_H4_max": float(d.alpha_H4.max()),
                "H_variants": json.dumps(
                    {k: int(v) for k, v in d.H.value_counts().sort_index().items()}
                ),
                "note": "coverage fraction over a declared deterministic envelope, not a probability"
                + (
                    "; EC is a stress test with no project source"
                    if env in ("EC", "EMC", "OAT-C")
                    else ""
                ),
            }
        )
    summary = pd.DataFrame(summ)
    summary.to_csv(OUT / "PCV05_summary.csv", index=False)
    write_json(
        OUT / "PCV05_meta.json",
        {
            "prereg_commit": "5d0b1986",
            "n_draws": len(df),
            "seed": SEED,
            "N_lhs": N_LHS,
            "machine_groups": {
                k: list(v[0]) + [v[1], v[2]] for k, v in MACHINE_GROUPS.items()
            },
            "converter_groups": {
                k: list(v[0]) + [v[1], v[2]] for k, v in CONV_GROUPS.items()
            },
            "nominal_alpha_H4": float(nominal.alpha_H4),
            "nominal_g_star": float(nominal.g_star),
            "elapsed_s": round(time.time() - started, 1),
        },
    )
    pd.set_option("display.width", 250)
    print(summary.drop(columns=["note"]).to_string(index=False))
    print("nominal alpha(H4) =", nominal.alpha_H4, "g* =", nominal.g_star)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
