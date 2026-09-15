"""Narrow TX4 TDS validation: one frozen D2 disturbance, three cases only.

This wrapper reuses the already validated G2 nonlinear phasor-domain solver;
it does not alter its disturbance, integration, observable, or outcome rules.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
OUT = ROOT / "results" / "TX4_TDS"
sys.path.insert(0, str(HERE.parent))
import G2_tds as g2  # noqa: E402


CASES = [
    ("IEEE-39", "proper 30+33+35", (30, 33, 35), {"g": 0.03625, "k": 1.425}, None),
    ("IEEE-39", "H4 original 30+33+35+37", (30, 33, 35, 37), {"g": 0.03625, "k": 1.425}, None),
    ("IEEE-39", "H4 retuned g=0.25", (30, 33, 35, 37), {"g": 0.25, "k": 1.425}, None),
]


def main() -> None:
    g2.pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.time()
    tasks = [(*case, "D2") for case in CASES]
    results = [g2.run(task) for task in tasks]
    rows = pd.DataFrame([r for r, _ in results])
    for row, trace in results:
        key = row["case"].replace(" ", "_").replace("+", "p").replace("=", "_")
        trace.to_csv(OUT / f"trace_{key}_D2.csv.gz", index=False)
    rows["f_err_hz"] = (rows.primary_f_hz - rows.lambda_f_hz).abs()
    rows["re_err"] = (rows.primary_re - rows.lambda_re).abs()
    rows["growth_sign_agrees"] = np.sign(rows.primary_re) == np.sign(rows.lambda_re)
    rows.to_csv(OUT / "TX4_TDS_FINAL.csv", index=False)
    summary = {
        "runs": int(len(rows)),
        "disturbance": g2.DISTURBANCES["IEEE-39"]["D2"],
        "pulse_percent": g2.STEP * 100.0,
        "pulse_duration_s": g2.PULSE_S,
        "horizon_rule": "clip(3/abs(lambda_re), 30, 120)",
        "observable_rule": "machine-speed pair maximizing critical-mode separation",
        "outcomes": {r["case"]: r["outcome"] for r, _ in results},
        "verdict_agreement": float(rows.verdict_agrees.mean()),
        "growth_sign_agreement": float(rows.growth_sign_agrees.mean()),
        "median_re_error": float(rows.re_err.median()),
        "max_re_error": float(rows.re_err.max()),
        "elapsed_s": time.time() - started,
        "model_scope": "nonlinear phasor-domain DAE; NOT EMT",
    }
    (OUT / "TX4_TDS_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(rows[["case", "lambda_re", "outcome", "primary_re", "lambda_f_hz", "primary_f_hz"]].to_string(index=False))


if __name__ == "__main__":
    main()
