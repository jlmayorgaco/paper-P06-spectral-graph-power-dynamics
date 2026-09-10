"""Regression: frozen F7 labels still reproduce with the post-freeze library.

Library files gained additive, default-off switches after the F7 freeze. This
re-solves random frozen F7 spot-check points on the direct path and compares
with the frozen direct labels.
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import evaluate
from _overnight import choose_workers, pin_blas_threads
from F7_policy_hypergraph import MAPS

N_PER_MAP = 40
SEED = 20260916


def _job(task):
    name, x, y, frozen = task
    return {
        "map": name,
        "frozen": frozen,
        "now": evaluate(MAPS[name].theta_xy(x, y), fast=False)["label"],
    }


def main() -> int:
    rng = np.random.default_rng(SEED)
    tasks = []
    for name in MAPS:
        s = pd.read_csv(RESULTS / "F7" / f"{name}_spotcheck.csv")
        for r in s.iloc[rng.choice(len(s), N_PER_MAP, replace=False)].itertuples():
            tasks.append((name, float(r.x), float(r.y), r.direct_label))
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        rows = pd.DataFrame(pool.map(_job, tasks, chunksize=2))
    out = {"points": int(len(rows)), "identical": int((rows.frozen == rows.now).sum())}
    (RESULTS / "F7_regression_after_freeze.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(out)
    return 0 if out["identical"] == out["points"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
