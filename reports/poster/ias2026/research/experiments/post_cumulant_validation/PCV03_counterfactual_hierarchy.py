# ruff: noqa: E501  -- table labels kept on one line
"""PCV03 - Phases 4 and 5: clean g-only counterfactual and the lower-order hierarchy.

Preregistration 5d0b1986, sections 6-7 and stopping rule S2. Reads PCV02 outputs for
the B0 truth and the baseline predictions at P4, G_S, G_S2, P_inf.

Counterfactual: P4 vs G_S (primary) and G_S2 (secondary) on k = 1.425, t = 1.5, h = 1;
only g changes. Identity checks: network Ybus, equilibrium voltages of H4, scheduled
P/Q at the candidate buses, ratings, static metrics. Direct bisection of g*.
Hierarchy: alpha space (B3, B4, B5, B5b, exact) and characteristic-function space
(local factors x Boolean truncations of det(I + Q) at orders 1..4; zero nearest the
exact critical eigenvalue by Newton).
"""

from __future__ import annotations

import hashlib
import sys

import numpy as np
import pandas as pd
from _pcv import (
    CORE,
    H4,
    POINTS,
    Realization,
    Static,
    label,
    out_dir,
    transverse_of,
    write_json,
)

from _f7_common import Theta, solve_subset
from ibr_cycles.cycles.connected import boolean_mobius, characteristic_values

OUT = out_dir("PCV03")
PCV02 = OUT.parent / "PCV02"
H_S = 1e-5
RADIUS = 1.0
FROZEN_GSTAR = 0.20768140450381395


def hash_array(a) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.round(np.asarray(a, complex), 10)).tobytes()
    ).hexdigest()[:16]


def alpha_h4(g):
    case = solve_subset(H4, Theta(g, 1.425, 1.5, 1.0))
    return transverse_of(case, classify=False)["alpha"]


def bisect_gstar(lo=0.2025, hi=0.2256, iters=40):
    a_lo, a_hi = alpha_h4(lo), alpha_h4(hi)
    assert a_lo > 0 > a_hi, (a_lo, a_hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if alpha_h4(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def truncated_char(r: Realization, s: complex, order: int) -> tuple[complex, complex]:
    """(local factor, Boolean truncation of det(I + Q_H4H4) at `order`) at s."""

    q, d, k, loc = r.q_matrix(s)
    vals = characteristic_values(q, [2] * 4)
    mu = boolean_mobius(vals, range(4))
    trunc = sum(v for t, v in mu.items() if len(t) <= order)
    return complex(np.linalg.det(loc)), complex(trunc)


def zero_of(r, lam0, order):
    def f(s):
        loc, tr = truncated_char(r, s, order)
        return loc * tr

    s = lam0
    for _ in range(60):
        v = f(s)
        dv = (f(s + H_S) - f(s - H_S)) / (2 * H_S)
        if abs(dv) < 1e-300:
            return None, "flat"
        step = v / dv
        s = s - step
        if abs(s - lam0) > RADIUS:
            return None, "no zero within radius"
        if abs(step) < 1e-12:
            return complex(s), "converged"
    return complex(s), "not converged"


def main(argv) -> int:
    truth = pd.read_csv(PCV02 / "PCV02_truth_core.csv")
    core = pd.read_csv(PCV02 / "PCV02_core_portfolios.csv")
    hyper = pd.read_csv(PCV02 / "PCV02_core_hypergraphs.csv")
    static = Static()

    # ---------------- S2: direct bisection of g* ---------------------------------------
    g_star = bisect_gstar()
    s2_ok = abs(g_star - FROZEN_GSTAR) <= 1e-4
    if not s2_ok:
        write_json(OUT / "PCV03_STOP.json", {"rule": "S2", "g_star": g_star})
        print("STOP S2", g_star)
        return 2

    # ---------------- identity checks across the counterfactual ------------------------
    cases = {
        pid: solve_subset(H4, Theta(*POINTS[pid]))
        for pid in ("P4", "G_S", "G_S2", "P_inf")
    }
    ident = {}
    net = cases["P4"].dae.network
    for pid, case in cases.items():
        v = case.dae.voltages(case.equilibrium.z)
        inj = {
            b: complex(
                case.dae.power_flow.injection(b, case.dae.network.ybus)
                + case.dae.network.loads.get(b, 0j)
            )
            for b in H4
        }
        ident[pid] = {
            "ybus_hash": hash_array(case.dae.network.ybus),
            "equilibrium_voltage_hash": hash_array(v),
            "max_voltage_diff_vs_P4": float(
                np.max(np.abs(v - cases["P4"].dae.voltages(cases["P4"].equilibrium.z)))
            ),
            "scheduled_P_MW": {b: round(inj[b].real * 100, 6) for b in H4},
            "scheduled_Q_Mvar": {b: round(inj[b].imag * 100, 6) for b in H4},
            "ratings_MVA": {b: float(net.machines[b]["Sn"]) for b in H4},
        }
    st = static.row(H4)

    rows = []

    def add(quantity, kind, values, note=""):
        same = len({str(values[p]) for p in ("P4", "G_S", "G_S2")}) == 1
        rows.append(
            {
                "quantity": quantity,
                "kind": kind,
                **{p: values[p] for p in ("P4", "G_S", "G_S2", "P_inf")},
                "identical_on_clean_line": same,
                "note": note,
            }
        )

    for q in ("g", "k", "t", "h"):
        i = "gkth".index(q)
        add(
            q,
            "policy coordinate",
            {p: POINTS[p][i] for p in POINTS},
            "only g changes on the clean line; P_inf also changes k",
        )
    add("network Ybus (hash)", "static", {p: ident[p]["ybus_hash"] for p in POINTS})
    add(
        "H4 equilibrium voltages (hash)",
        "static",
        {p: ident[p]["equilibrium_voltage_hash"] for p in POINTS},
        "max abs difference vs P4: "
        + ", ".join(f"{p} {ident[p]['max_voltage_diff_vs_P4']:.1e}" for p in POINTS),
    )
    add(
        "scheduled P at 30/33/35/37 [MW]",
        "static",
        {p: str(ident[p]["scheduled_P_MW"]) for p in POINTS},
    )
    add(
        "scheduled Q at 30/33/35/37 [Mvar]",
        "static",
        {p: str(ident[p]["scheduled_Q_Mvar"]) for p in POINTS},
    )
    add("ratings [MVA]", "static", {p: str(ident[p]["ratings_MVA"]) for p in POINTS})
    for key in (
        "pg_mw",
        "sn_mva",
        "penetration",
        "gscr",
        "min_scr",
        "max_miif",
        "compactness_score",
    ):
        add(
            f"B1/B2/B7 {key}",
            "static (controller-independent by construction)",
            {p: round(st[key], 10) for p in POINTS},
        )
    dyn = {}
    for pid in POINTS:
        d = core[(core.point == pid) & (core.subset == label(H4))].iloc[0]
        h = hyper[hyper.point == pid].iloc[0]
        tr = truth[(truth.point == pid) & (truth.subset == label(H4))].iloc[0]
        dyn[pid] = {
            "alpha": d.alpha,
            "crit_hz": d.crit_hz,
            "rhp": tr.rhp,
            "H": h.H_true,
            "kappa": h.kappa_true,
            "status": d.status,
            "B3": d.B3_alpha,
            "B4": d.B4_alpha,
            "B5": d.B5_alpha,
            "B5b": d.B5b_alpha,
            "closure_nonoracle": d.closure_nonoracle,
            "closure_oracle": d.closure_oracle,
            "abs_chi_nonoracle": d.abs_chi_nonoracle,
        }
    for key, kind in (
        ("alpha", "dynamic (B0)"),
        ("crit_hz", "dynamic (B0)"),
        ("rhp", "dynamic (B0)"),
        ("status", "dynamic (B0)"),
        ("H", "dynamic (B0)"),
        ("kappa", "dynamic (B0)"),
        ("B3", "baseline B3 alpha-hat(H4)"),
        ("B4", "baseline B4 alpha-hat(H4)"),
        ("B5", "baseline B5 alpha-hat(H4)"),
        ("B5b", "baseline B5b alpha-hat(H4)"),
        ("closure_nonoracle", "B8a"),
        ("closure_oracle", "B8b"),
        ("abs_chi_nonoracle", "B9a (secondary)"),
    ):
        add(f"{key}", kind, {p: dyn[p][key] for p in POINTS})
    add(
        "g* on the clean line (direct bisection)",
        "boundary",
        {p: g_star for p in POINTS},
        f"frozen {FROZEN_GSTAR}",
    )
    cf = pd.DataFrame(rows)
    cf.to_csv(OUT / "PCV03_counterfactual.csv", index=False)

    # ---------------- hierarchy ----------------------------------------------------------
    hier = []
    for pid in ("P4", "G_S"):
        th = POINTS[pid]
        tr = truth[(truth.point == pid) & (truth.subset == label(H4))].iloc[0]
        lam = complex(tr.crit_re, 2 * np.pi * tr.crit_hz)
        d = core[(core.point == pid) & (core.subset == label(H4))].iloc[0]
        for level, val in (
            ("B3 first-order modal sensitivity", d.B3_alpha),
            ("B4 additive singles", d.B4_alpha),
            ("B5 pairwise", d.B5_alpha),
            ("B5b third order", d.B5b_alpha),
            ("exact (B0)", d.alpha),
        ):
            hier.append(
                {
                    "point": pid,
                    "space": "alpha",
                    "level": level,
                    "predicted_alpha": val,
                    "predicted_unstable": bool(val > 0),
                    "residual_to_exact": val - d.alpha,
                    "flags_H4": bool(val > 0),
                }
            )
        r = Realization(CORE, th)
        for order, name in (
            (1, "local factors only (order <= 1)"),
            (2, "pairwise closure (order <= 2)"),
            (3, "third-order closure (order <= 3)"),
            (4, "exact full closure (order 4)"),
        ):
            loc, trv = truncated_char(r, lam, order)
            z, status = zero_of(r, lam, order)
            hier.append(
                {
                    "point": pid,
                    "space": "characteristic",
                    "level": name,
                    "predicted_alpha": None if z is None else z.real,
                    "predicted_freq_hz": None if z is None else z.imag / (2 * np.pi),
                    "zero_status": status,
                    "predicted_unstable": None if z is None else bool(z.real > 0),
                    "abs_truncation_at_exact_lambda": abs(trv),
                    "abs_local_at_exact_lambda": abs(loc),
                    "residual_to_exact": None if z is None else abs(z - lam),
                    "flags_H4": None if z is None else bool(z.real > 0),
                }
            )
        hier.append(
            {
                "point": pid,
                "space": "reference",
                "level": "exact critical eigenvalue",
                "predicted_alpha": lam.real,
                "predicted_freq_hz": lam.imag / (2 * np.pi),
            }
        )
    hier_df = pd.DataFrame(hier)
    hier_df.to_csv(OUT / "PCV03_hierarchy.csv", index=False)
    write_json(
        OUT / "PCV03_summary.json",
        {
            "g_star_direct": g_star,
            "g_star_frozen": FROZEN_GSTAR,
            "S2_ok": s2_ok,
            "identity": ident,
            "static_H4": st,
            "dynamic_H4": dyn,
        },
    )
    pd.set_option("display.width", 250)
    print(cf.to_string(index=False))
    print(hier_df.to_string(index=False))
    print("g* =", g_star)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
