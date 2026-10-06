"""Reduced-only g-boundary search for the TX4 H4 portfolio.

No full-order eigenvalues or historical answer tables are imported here.  A
fresh baseline/local port library is built for each g evaluation and the
boundary is located by minimizing the collective closure singular value over
the frozen frequency band.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
EXP = HERE.parent
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(EXP.parent / "src"))
from _f7_common import Theta  # noqa: E402
import tx4_blind_predict as bp  # noqa: E402

OUT = ROOT / "results"
CORE = (30, 33, 35, 37)
FREQ_GRID = np.linspace(0.3, 1.5, 81)


def evaluate(g: float):
    bp.THETA = Theta(g=float(g), k=1.425, t=1.5, h=1.0)
    kernel = bp.make_kernel(CORE)
    values = np.array([
        kernel.qdata(CORE, 2j * np.pi * float(freq))["collective_sigma_min"]
        for freq in FREQ_GRID
    ])
    j = int(np.argmin(values))
    lo = float(FREQ_GRID[max(0, j - 1)])
    hi = float(FREQ_GRID[min(len(FREQ_GRID) - 1, j + 1)])
    if hi > lo:
        opt = minimize_scalar(
            lambda f: kernel.qdata(CORE, 2j * np.pi * float(f))["collective_sigma_min"],
            bounds=(lo, hi), method="bounded", options={"xatol": 1e-10},
        )
        freq, sigma = float(opt.x), float(opt.fun)
    else:
        freq, sigma = float(FREQ_GRID[j]), float(values[j])
    q = kernel.qdata(CORE, 2j * np.pi * freq)
    return {"g": float(g), "frequency_hz": freq, "collective_sigma_min": sigma,
            "local_sigma_min": q["local_sigma_min"],
            "local_det_min_abs": q["local_det_min_abs"]}


def main():
    result = minimize_scalar(
        lambda g: evaluate(float(g))["collective_sigma_min"],
        bounds=(0.19, 0.225), method="bounded",
        options={"xatol": 1e-9, "maxiter": 40},
    )
    center = evaluate(float(result.x))
    side_lo = evaluate(float(result.x - 1e-4))
    side_hi = evaluate(float(result.x + 1e-4))
    rows = pd.DataFrame([side_lo, center, side_hi])
    rows["search"] = "reduced_port_collective_sigma_min"
    rows["protocol_status"] = "REDUCED_REPRODUCTION_AFTER_REVEAL"
    rows.to_csv(OUT / "TX4_G_BOUNDARY_BLIND_PREDICTION.csv", index=False)
    summary = {
        "g_return_reduced": center["g"],
        "frequency_return_hz": center["frequency_hz"],
        "collective_sigma_min": center["collective_sigma_min"],
        "local_sigma_min": center["local_sigma_min"],
        "violating_side_checks": [side_lo, side_hi],
        "protocol_status": "REDUCED_REPRODUCTION_AFTER_REVEAL",
        "answer_key_loaded_by_this_script": False,
    }
    (OUT / "TX4_G_BOUNDARY_REDUCED_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
