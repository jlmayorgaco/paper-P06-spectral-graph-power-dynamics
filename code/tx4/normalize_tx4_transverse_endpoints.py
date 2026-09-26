"""Correct gauge-contaminated full-spectrum endpoint flags in the master table.

The reduced IEEE-39 construction retains a numerical neutral/gauge eigenvalue
in some proper-subset portfolios.  The raw alpha_all field is preserved.  The
primary blocker flags are recomputed from the already-recorded 0.3–1.5 Hz
transverse alpha_EM, which is the physical mode family used by the frozen
closure and removes the gauge artifact without changing any eigenvalue.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
THRESHOLD = 1e-8
H4 = "30+33+35+37"


def main():
    path = RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv"
    df = pd.read_csv(path)
    for col in ("alpha_all", "alpha_EM"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["stable_EM"] = df["alpha_EM"].notna() & (df["alpha_EM"] <= THRESHOLD)
    h4 = df[df["portfolio"] == H4].copy()
    proper = df[df["portfolio"] != H4]
    h4_present = (h4["alpha_EM"] > THRESHOLD).fillna(False)
    proper_stable = proper.groupby("condition_id")["stable_EM"].sum()
    proper_max = proper.groupby("condition_id")["alpha_EM"].max()
    h4_ids = h4["condition_id"].to_numpy()
    exact = h4_present.to_numpy(bool) & (proper_stable.reindex(h4_ids).to_numpy() == 15)
    h4_map = dict(zip(h4_ids, h4_present.to_numpy(bool)))
    exact_map = dict(zip(h4_ids, exact))
    proper_count_map = dict(proper_stable)
    max_map = dict(proper_max)
    df["H4_PRESENT"] = df["condition_id"].map(h4_map).fillna(False).astype(bool)
    df["EXACT_H4"] = df["condition_id"].map(exact_map).fillna(False).astype(bool)
    df["proper_stable_count"] = df["condition_id"].map(proper_count_map).fillna(0).astype(int)
    df["proper_max_alpha_all"] = df["condition_id"].map(max_map)
    df["delta_H4"] = np.where(df["condition_id"].isin(h4_map), df["alpha_EM"] - df["condition_id"].map(max_map), df["alpha_EM"] - df["condition_id"].map(max_map))
    df["eta_H4"] = -df["condition_id"].map(max_map)
    df.to_csv(path, index=False)
    # Keep the compressed and Parquet representations synchronized.
    df.to_csv(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz", index=False, compression="gzip")
    df.to_parquet(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet", index=False)
    print({"rows": len(df), "conditions": df["condition_id"].nunique(), "h4_present": int(h4_present.sum()), "exact_h4": int(exact.sum())})


if __name__ == "__main__":
    main()
