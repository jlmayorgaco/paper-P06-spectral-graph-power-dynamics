# ruff: noqa: E501  -- case tables kept on one line
"""EMT-PRED - frozen PHASOR predictions for every preregistered EMT case (tx3-analysis).

Run BEFORE any TX4 EMT result; committed with the EMT preregistration
(docs/20260911_PAREMT_EMT_PREREG_V1.md). Writes
  results/EMT_PRED/emt_cases.json        the case definitions the EMT runner consumes
  results/EMT_PRED/phasor_predictions.csv the frozen phasor prediction of every case
  results/EMT_PRED/emt16_selection.json   the deterministic robustness-holdout draw selection

Prediction per case (canonical phasor model, direct path):
  alpha_perp  maximum real part of the transverse spectrum (frozen definition)
  status      FC01 / FC03 logic (sign count if min|Re| >= 2e-3, else the h/2h classifier)
  band_re, band_hz  rightmost transverse eigenvalue with 0.2 <= f <= 1.2 Hz (the EMT estimator band)
Governed cases use the governed quotient C = span{R_x} (FC03). Condenser cases use G1 R_0010000.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
for p in (
    EXP,
    EXP / "final_closure",
    EXP / "post_cumulant_validation",
    EXP / "connected_cumulants",
):
    sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _bootstrap  # noqa: E402,F401
from _f7_common import Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from G1_f8_rhp_followup import _em_config  # noqa: E402
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.governed import govern  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402

RESEARCH = HERE.parents[1]
OUT = RESEARCH / "results" / "EMT_PRED"
H4 = (30, 33, 35, 37)
P4 = (0.03625, 1.425, 1.5, 1.0)
G_S = (0.25, 1.425, 1.5, 1.0)
P_INF = (1.0, 0.5, 1.5, 1.0)
BAND = (0.2, 1.2)
MARGIN = 2e-3
NEWTON_G1 = 0.08806432283854684  # frozen FC18 first port-Newton iterate from P4
G_STAR_PHASOR = 0.20768140450381395
SEL_SEED = 20260921
PCV05 = RESEARCH / "results" / "PCV" / "PCV05" / "PCV05_draws.csv"
HOLDOUT = (
    RESEARCH
    / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff/dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation/holdout_lines.json"
)
_CFG: dict = {}


def label(m):
    return "+".join(map(str, sorted(m))) or "BASE"


# ------------------------------------------------------------------ cases --
def cases() -> list[dict]:
    out = []

    def add(cid, exp, members, theta, variant=None, note="", **extra):
        out.append(
            {
                "case_id": cid,
                "experiment": exp,
                "members": list(members),
                "theta": list(theta),
                "variant": variant or {"kind": "none"},
                "note": note,
                **extra,
            }
        )

    subsets = [tuple(s) for r in range(5) for s in combinations(H4, r)]
    for s in subsets:
        add(f"EMT05_P4_{label(s)}", "EMT05", s, P4)
    add(
        "EMT04_G_S_H4",
        "EMT04",
        H4,
        G_S,
        note="stable near-policy point (dt convergence)",
    )
    for g in (0.18, 0.20, 0.205, 0.21, 0.225, 0.25):
        add(f"EMT08_g{g:.3f}", "EMT08", H4, (g, 1.425, 1.5, 1.0))
    add("EMT09_P_inf_H4", "EMT09", H4, P_INF, note="multi-coordinate (g and k change)")
    add(
        "EMT10_newton1",
        "EMT10",
        H4,
        (NEWTON_G1, 1.425, 1.5, 1.0),
        note="frozen first port-Newton iterate (under-correction)",
    )
    add(
        "EMT10_gstar",
        "EMT10",
        H4,
        (G_STAR_PHASOR, 1.425, 1.5, 1.0),
        note="phasor Newton boundary (near zero; UNRESOLVED allowed)",
    )
    for k in (1.25, 1.28, 1.30, 1.31, 1.33, 1.35):
        add(f"EMT11_k{k:.2f}", "EMT11", H4, (0.03625, k, 1.5, 1.0))
    lines = json.loads(HOLDOUT.read_text())["lines"]
    for li in lines:
        add(
            f"EMT12_L{int(li):02d}",
            "EMT12",
            H4,
            P4,
            {"kind": "line", "line_index": int(li), "gamma": 1.5},
        )
    for p in (0.020, 0.025, 0.030, 0.050):
        add(
            f"EMT13_cond{100 * p:.1f}pct",
            "EMT13",
            H4,
            P4,
            {"kind": "condenser", "config": "R_0010000", "fraction": p},
        )
    add("EMT14_gov_P4_H4", "EMT14", H4, P4, {"kind": "governed"})
    add(
        "EMT14_gov_g0.020_H4",
        "EMT14",
        H4,
        (0.020, 1.425, 1.5, 1.0),
        {"kind": "governed"},
        note="inside the frozen governed kappa=4 interval g in [0.016, 0.023]",
    )
    for g in (0.22, 0.30):
        add(
            f"EMT15_g{g:.2f}",
            "EMT15",
            H4,
            (g, 1.425, 1.5, 1.0),
            note="amplitude sweep point",
        )
    sel = emt16_selection()
    for d in sel["draws"]:
        wit = tuple(d["witness"]) if d["witness"] else H4
        tests = (
            [wit] + [tuple(c) for c in combinations(wit, len(wit) - 1)]
            if d["group"] in ("A", "B")
            else [H4]
        )
        tests.append(())
        for m in tests:
            add(
                f"EMT16_{d['group']}{d['rank']}_{d['envelope']}_{d['draw']}_{label(m)}",
                "EMT16",
                m,
                P4,
                {
                    "kind": "draw",
                    "envelope": d["envelope"],
                    "draw": d["draw"],
                    "factors": d["factors"],
                },
            )
    return out


def emt16_selection() -> dict:
    df = pd.read_csv(PCV05)
    df = df[df.exact.astype(bool)]
    mach = df[df.envelope.isin(["EM-f", "EM-u", "EMC"])]
    groups = {
        "A": mach[mach.H == "30+33+35+37"],
        "B": mach[mach.H == "30+33+35"],
        "C": mach[mach.H == "EMPTY"],
    }
    n = {"A": 6, "B": 6, "C": 4}
    rng = np.random.default_rng(SEL_SEED)
    draws = []
    for g in ("A", "B", "C"):
        pool = groups[g].sort_values(["envelope", "draw"]).reset_index(drop=True)
        idx = sorted(rng.choice(len(pool), size=n[g], replace=False).tolist())
        for rank, i in enumerate(idx):
            r = pool.iloc[i]
            wit = [int(b) for b in r.H.split("+")] if r.H != "EMPTY" else []
            draws.append(
                {
                    "group": g,
                    "rank": rank,
                    "envelope": r.envelope,
                    "draw": int(r.draw),
                    "H": r.H,
                    "witness": wit,
                    "factors": json.loads(r.factors),
                }
            )
    return {
        "seed": SEL_SEED,
        "rule": "A: H={H4}; B: H={30+33+35} (single hyperedge); C: H=EMPTY; "
        "machine envelopes EM-f/EM-u/EMC, exact draws only, sorted by (envelope, draw), rng.choice without replacement",
        "eligible_counts": {g: int(len(v)) for g, v in groups.items()},
        "draws": draws,
    }


# ---------------------------------------------------------- predictions --
def _reduced(j):
    return j.fx - j.fz @ np.linalg.solve(j.gz, j.gx)


def spectrum_status(dae, x, z, governed):
    j1 = central_difference_jacobians(dae, x, z, {})
    j2 = central_difference_jacobians(dae, x, z, {}, scale_x=2.0, scale_z=2.0)
    a = _reduced(j1)
    d = a - _reduced(j2)
    r_x, _ = rotation_generator(dae, z)
    w = None if governed else frequency_partner(dae).w
    dead = np.flatnonzero(~np.any(a != 0.0, axis=1))
    if dead.size:
        keep = np.setdiff1d(np.arange(a.shape[0]), dead)
        a, d, r_x = a[np.ix_(keep, keep)], d[np.ix_(keep, keep)], r_x[keep]
        w = None if w is None else w[keep]
    tr = transverse_operator(a, r_x, w)
    ev = np.linalg.eigvals(tr.a_perp)
    if np.abs(ev.real).min() < MARGIN:
        st = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    else:
        st = "UNSTABLE" if (ev.real > 0).any() else "STABLE"
    return st, ev


def scaled_network(li, gamma):
    payload = json.loads((RESEARCH / "configs/ias2026/ieee39_network.json").read_text())
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {b: i for i, b in enumerate(order)}
    y = np.zeros((len(order), len(order)), complex)
    for ell, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        s = gamma if ell == li else 1.0
        series = s / complex(line["r"], line["x"])
        charging = s * complex(line["g"], line["b"]) / 2
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        y[f, f] += (series + charging) / abs(m) ** 2
        y[t, t] += series + charging
        y[f, t] += -series / np.conj(m)
        y[t, f] += -series / m
    for sh in payload["shunts"]:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return replace(
        load_network(RESEARCH / "configs/ias2026/ieee39_network.json"), ybus=y
    )


def draw_kwargs(factors, theta):
    sys.path.insert(0, str(EXP / "post_cumulant_validation"))
    from PCV05_robustness import case_kwargs

    draw = {
        "fleet": factors["fleet"],
        "unit": {int(b): v for b, v in factors["unit"].items()},
        "conv": {int(b): v for b, v in factors["conv"].items()},
    }
    return case_kwargs(draw, theta)


def predict(case):
    pin_blas_threads()
    if not _CFG:
        _CFG.update({c["name"]: c for c in r_configs()})
    v = case["variant"]
    th = Theta(*case["theta"])
    members = tuple(case["members"])
    governed = False
    if v["kind"] == "condenser":
        c = solve_config(members, th, _em_config(_CFG[v["config"]], v["fraction"]))
        dae, x, z = c.dae, c.equilibrium.x, c.equilibrium.z
    elif v["kind"] == "governed":
        c = solve_config(members, th, _CFG["R_none"])
        dae, x, z = govern(c)
        governed = True
    elif v["kind"] == "line":
        from _cc import kwargs_of

        c = solve_case(
            ReplacementPlan.of({b: 1.0 for b in members}),
            network=scaled_network(v["line_index"], v["gamma"]),
            **kwargs_of(tuple(case["theta"])),
        )
        dae, x, z = c.dae, c.equilibrium.x, c.equilibrium.z
    elif v["kind"] == "draw":
        c = solve_case(
            ReplacementPlan.of({b: 1.0 for b in members}),
            **draw_kwargs(v["factors"], tuple(case["theta"])),
        )
        dae, x, z = c.dae, c.equilibrium.x, c.equilibrium.z
    else:
        c = solve_config(members, th, _CFG["R_none"])
        dae, x, z = c.dae, c.equilibrium.x, c.equilibrium.z
    st, ev = spectrum_status(dae, x, z, governed)
    f = ev.imag / (2 * np.pi)
    band = ev[(f >= BAND[0]) & (f <= BAND[1])]
    crit = complex(band[np.argmax(band.real)]) if band.size else complex("nan")
    return {
        "case_id": case["case_id"],
        "experiment": case["experiment"],
        "members": label(members),
        "g": case["theta"][0],
        "k": case["theta"][1],
        "t": case["theta"][2],
        "h": case["theta"][3],
        "variant": json.dumps({k: val for k, val in v.items() if k != "factors"}),
        "status": st,
        "alpha_perp": float(ev.real.max()),
        "n_rhp": int((ev.real > 0).sum()),
        "min_abs_re": float(np.abs(ev.real).min()),
        "band_re": crit.real,
        "band_hz": crit.imag / (2 * np.pi),
        "note": case["note"],
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cs = cases()
    (OUT / "emt_cases.json").write_text(
        json.dumps({"band_hz": BAND, "cases": cs}, indent=1), encoding="utf-8"
    )
    (OUT / "emt16_selection.json").write_text(
        json.dumps(emt16_selection(), indent=1), encoding="utf-8"
    )
    with Pool(16) as pool:
        rows = pool.map(predict, cs, chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "phasor_predictions.csv", index=False)
    pd.set_option("display.width", 250)
    print(
        df[
            ["case_id", "status", "alpha_perp", "band_re", "band_hz", "n_rhp"]
        ].to_string(index=False)
    )
    print(len(df), "cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
