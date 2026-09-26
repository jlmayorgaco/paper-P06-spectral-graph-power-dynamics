# ruff: noqa: E501  -- table labels kept on one line
"""PCV06 step 3: compare the ANDES reproduction with the internal reference.

Spec: docs/20260911_GFL_REPRODUCTION_SPEC.md (commit 3e847a4f), sections 7-9.
Reads the ANDES outputs and the internal reference (JSON only; imports nothing from
either tool).

    python compare.py gate0   -> PCV06_gate0.csv, PCV06_gate0_summary.json
    python compare.py cases   -> PCV06_cases.csv, PCV06_comparison.json,
                                 results/20260911_GFL_REPRODUCTION.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[2]
OUT = RESEARCH / "results" / "PCV" / "PCV06"
GATE0_TOL = 1e-9
TOL = {
    "voltage": 1e-6,
    "residual": 1e-6,
    "alpha": 1e-3,
    "freq": 1e-3,
    "band_nearest": 1e-4,
}
BAND = (0.3, 1.5)


def load(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def gate0() -> int:
    ref = load("PCV06_internal_reference.json")["gate0_outputs"]
    got = load("PCV06_andes_gate0.json")
    rows = []
    for r, a in zip(ref, got["rows"], strict=True):
        assert (r["kind"], r["bus"], r["k"]) == (a["kind"], a["bus"], a["k"])
        fi, fa = np.array(r["f"]), np.array(a["f"])
        scale = max(1.0, float(np.abs(fi).max()))
        rows.append(
            {
                "kind": r["kind"],
                "bus": r["bus"],
                "k": r["k"],
                "f_rel_err": float(np.abs(fa - fi).max() / scale),
                "worst_state": int(np.argmax(np.abs(fa - fi))),
                "P_abs_err": abs(a["P_inj"] - r["P_inj"]) / max(1.0, abs(r["P_inj"])),
                "Q_abs_err": abs(a["Q_inj"] - r["Q_inj"]) / max(1.0, abs(r["Q_inj"])),
            }
        )
    df = pd.DataFrame(rows)
    df["pass"] = (
        (df.f_rel_err <= GATE0_TOL)
        & (df.P_abs_err <= GATE0_TOL)
        & (df.Q_abs_err <= GATE0_TOL)
    )
    df.to_csv(OUT / "PCV06_gate0.csv", index=False)
    per = (
        df.groupby(["kind", "bus"])
        .agg(
            f=("f_rel_err", "max"),
            P=("P_abs_err", "max"),
            Q=("Q_abs_err", "max"),
            ok=("pass", "all"),
        )
        .reset_index()
    )
    summary = {
        "tolerance": GATE0_TOL,
        "n_vectors": len(df),
        "all_pass": bool(df["pass"].all()),
        "max_f_rel_err": float(df.f_rel_err.max()),
        "max_P_err": float(df.P_abs_err.max()),
        "max_Q_err": float(df.Q_abs_err.max()),
        "per_instance": per.to_dict(orient="records"),
        "andes_init_residual_base": got["residual_base"],
        "andes_init_residual_h4": got["residual_h4"],
    }
    (OUT / "PCV06_gate0_summary.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )
    print(per.to_string(index=False))
    print(
        "GATE 0",
        "PASS" if summary["all_pass"] else "FAIL",
        summary["max_f_rel_err"],
        summary["max_P_err"],
        summary["max_Q_err"],
    )
    return 0 if summary["all_pass"] else 3


def hypergraph(verdicts: dict) -> tuple[str, float]:
    if verdicts["BASE"] != "STABLE":
        return f"BASE_{verdicts['BASE']}", float("nan")
    unsafe = [
        tuple(int(b) for b in s.split("+"))
        for s, v in verdicts.items()
        if v == "UNSTABLE"
    ]
    minimal = sorted(
        (s for s in unsafe if not any(set(r) < set(s) for r in unsafe)),
        key=lambda s: (len(s), s),
    )
    return "|".join("+".join(map(str, e)) for e in minimal) or "EMPTY", float(
        min((len(e) for e in minimal), default=np.inf)
    )


def cases() -> int:
    ref = {
        (c["point"], c["subset"]): c
        for c in load("PCV06_internal_reference.json")["cases"]
    }
    got = load("PCV06_andes_cases.json")
    rows = []
    for a in got["cases"]:
        r = ref[(a["point"], a["subset"])]
        ev_i = np.array([complex(*z) for z in r["spectrum"]])
        ev_a = np.array([complex(*z) for z in a["spectrum"]])
        f_i = ev_i.imag / (2 * np.pi)
        band_i = ev_i[(f_i >= BAND[0]) & (f_i <= BAND[1])]
        nearest = max((float(np.abs(ev_a - z).min()) for z in band_i), default=0.0)
        dv = max(
            abs(complex(*r["voltages"][b]) - complex(*a["voltages"][b]))
            for b in r["voltages"]
        )
        rows.append(
            {
                "point": a["point"],
                "subset": a["subset"],
                "internal_status": r["status"],
                "andes_verdict": a["verdict"],
                "internal_alpha": r["alpha"],
                "andes_alpha": a["alpha"],
                "delta_alpha": a["alpha"] - r["alpha"],
                "internal_crit_hz": r["crit_hz"],
                "andes_crit_hz": a["crit_hz"],
                "delta_freq_hz": a["crit_hz"] - r["crit_hz"],
                "internal_rhp": r["rhp"],
                "andes_rhp": a["rhp"],
                "dim_internal_perp": r["n_x"] - 2,
                "dim_andes_transverse": len(ev_a),
                "andes_init_residual": a["init_residual"],
                "max_voltage_diff": dv,
                "structural_pair_ok": a["structural_pair_ok"],
                "structural_pair_moduli": json.dumps(a["structural_pair"]),
                "band_nearest_max": nearest,
                "andes_min_abs_re": a["min_abs_re"],
            }
        )
    df = pd.DataFrame(rows)
    df["pass_alpha"] = df.delta_alpha.abs() <= TOL["alpha"]
    df["pass_freq"] = df.delta_freq_hz.abs() <= TOL["freq"]
    df["pass_rhp"] = df.internal_rhp == df.andes_rhp
    df["pass_verdict"] = df.internal_status == df.andes_verdict
    df["pass_equilibrium"] = (df.andes_init_residual <= TOL["residual"]) & (
        df.max_voltage_diff <= TOL["voltage"]
    )
    df["pass_dimension"] = df.dim_internal_perp == df.dim_andes_transverse
    df["pass_band_nearest_secondary"] = df.band_nearest_max <= TOL["band_nearest"]
    df.to_csv(OUT / "PCV06_cases.csv", index=False)
    hyper = {}
    for pid in df.point.unique():
        d = df[df.point == pid]
        hi = hypergraph(dict(zip(d.subset, d.internal_status, strict=True)))
        ha = hypergraph(dict(zip(d.subset, d.andes_verdict, strict=True)))
        hyper[pid] = {
            "internal_H": hi[0],
            "internal_kappa": hi[1],
            "andes_H": ha[0],
            "andes_kappa": ha[1],
            "andes_script_H": got["hypergraphs"][pid]["H"],
            "equal": hi == ha,
        }
    primary = {
        k: bool(df[k].all())
        for k in (
            "pass_alpha",
            "pass_freq",
            "pass_rhp",
            "pass_verdict",
            "pass_equilibrium",
        )
    }
    primary["pass_H_kappa"] = all(h["equal"] for h in hyper.values())
    stop = bool(
        (~df.pass_verdict).any()
        or (df.delta_alpha.abs() > 1e-2).any()
        or not primary["pass_H_kappa"]
    )
    summary = {
        "tolerances": TOL,
        "n_cases": len(df),
        "primary": primary,
        "gate5_pass": all(primary.values()),
        "stop_rule_triggered": stop,
        "max_abs_delta_alpha": float(df.delta_alpha.abs().max()),
        "max_abs_delta_freq_hz": float(df.delta_freq_hz.abs().max()),
        "max_voltage_diff": float(df.max_voltage_diff.max()),
        "max_init_residual": float(df.andes_init_residual.max()),
        "structural_pair_ok_all": bool(df.structural_pair_ok.all()),
        "dimension_ok_all": bool(df.pass_dimension.all()),
        "secondary_band_nearest_ok": bool(df.pass_band_nearest_secondary.all()),
        "max_band_nearest": float(df.band_nearest_max.max()),
        "hypergraphs": hyper,
    }
    (OUT / "PCV06_comparison.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )
    cols = [
        "point",
        "subset",
        "internal_status",
        "andes_verdict",
        "internal_alpha",
        "andes_alpha",
        "delta_alpha",
        "internal_crit_hz",
        "andes_crit_hz",
        "delta_freq_hz",
        "internal_rhp",
        "andes_rhp",
        "max_voltage_diff",
        "andes_init_residual",
        "band_nearest_max",
        "pass_alpha",
        "pass_freq",
        "pass_rhp",
        "pass_verdict",
        "pass_equilibrium",
        "pass_band_nearest_secondary",
    ]
    df[cols].to_csv(RESEARCH / "results" / "20260911_GFL_REPRODUCTION.csv", index=False)
    pd.set_option("display.width", 250)
    print(
        df[
            [
                "point",
                "subset",
                "internal_status",
                "andes_verdict",
                "delta_alpha",
                "delta_freq_hz",
                "max_voltage_diff",
                "andes_init_residual",
                "band_nearest_max",
            ]
        ].to_string(index=False)
    )
    print(
        json.dumps({k: v for k, v in summary.items() if k != "hypergraphs"}, indent=1)
    )
    print(json.dumps(hyper, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(gate0() if sys.argv[1] == "gate0" else cases())
