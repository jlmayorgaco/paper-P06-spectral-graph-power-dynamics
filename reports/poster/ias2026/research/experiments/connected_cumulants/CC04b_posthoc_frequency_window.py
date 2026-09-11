# ruff: noqa: E501  -- report line kept on one line
"""CC04b - POST HOC, DESCRIPTIVE ONLY (not preregistered; no verdict depends on it).

After CC04 it was noticed that the failing portfolios' critical band modes lie at
0.42-0.57 Hz while most stable controls' lie at 0.7-1.4 Hz. This script repeats the
CC04 AUCs (same orientations) inside the 0.35-0.60 Hz window only, to see whether the
separation is a frequency artifact. It reads results/CC/CC04/CC04_portfolios.csv.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parents[2] / "results" / "CC" / "CC04"
ORIENT = {
    "abs_chi_axis": 1,
    "abs_mu_axis": 1,
    "cycle_score_axis": 1,
    "closure_distance_axis": -1,
    "nu_at_mode": 1,
    "pg_mw": 1,
    "sn_mva": 1,
    "gscr": -1,
}


def auc(x, labels):
    p, n = x[labels], x[~labels]
    return float(
        ((p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum())
        / (len(p) * len(n))
    )


def main() -> int:
    df = pd.read_csv(HERE / "CC04_portfolios.csv")
    w = df[(df.lambda_hz >= 0.35) & (df.lambda_hz <= 0.60)]
    labels = (w.role == "failing").to_numpy()
    lines = [
        f"POST HOC, DESCRIPTIVE. window 0.35-0.60 Hz: failing {labels.sum()}, controls {(~labels).sum()}"
    ]
    for name, sign in ORIENT.items():
        lines.append(f"{name} {auc(sign * w[name].to_numpy(float), labels):.3f}")
    rho = spearmanr(df.abs_chi_axis, df.closure_distance_axis).statistic
    lines.append(f"spearman |chi| vs closure distance, all 35: {rho:.3f}")
    text = "\n".join(lines)
    (HERE / "CC04_posthoc_frequency_window.txt").write_text(
        text + "\n", encoding="utf-8"
    )
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
