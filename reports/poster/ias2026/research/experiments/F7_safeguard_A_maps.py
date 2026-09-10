"""F7 Safeguard A over the maps: structure checks at the spot-check points, direct path.

At every spot-check point of every map and for every subset: the state dimension,
the number of eigenvalues below 1e-3 (the structural angle-reference pair only,
so exactly two), and cond(g_z).
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import LABELS, SUBSETS, solve_subset
from _overnight import choose_workers, pin_blas_threads
from F7_policy_hypergraph import MAPS

F7 = RESULTS / "F7"


def _check(theta):
    out = {}
    for members, label in zip(SUBSETS, LABELS, strict=True):
        case = solve_subset(members, theta)
        values = np.linalg.eigvals(case.system.A)
        out[f"nx_{label}"] = case.n_states
        out[f"small_{label}"] = int(np.sum(np.abs(values) < 1e-3))
        out[f"gz_{label}"] = case.gz_condition
    return out


def main() -> int:
    report = {}
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        for name, spec in MAPS.items():
            spots = pd.read_csv(F7 / f"{name}_spotcheck.csv")
            thetas = [
                spec.theta_xy(float(x), float(y))
                for x, y in zip(spots.x, spots.y, strict=True)
            ]
            rows = pd.DataFrame(pool.map(_check, thetas, chunksize=2))
            report[name] = {
                "points": len(rows),
                "state_dimension_constant_per_subset": bool(
                    all(rows[f"nx_{lab}"].nunique() == 1 for lab in LABELS)
                ),
                "state_dimension_by_subset": {
                    lab: int(rows[f"nx_{lab}"].iloc[0]) for lab in LABELS
                },
                "small_eigenvalue_counts_seen": sorted(
                    set(int(v) for lab in LABELS for v in rows[f"small_{lab}"].unique())
                ),
                "gz_condition_max_relative_spread": float(
                    max(
                        np.ptp(rows[f"gz_{lab}"]) / rows[f"gz_{lab}"].median()
                        for lab in LABELS
                    )
                ),
                "gz_condition_range": [
                    float(min(rows[f"gz_{lab}"].min() for lab in LABELS)),
                    float(max(rows[f"gz_{lab}"].max() for lab in LABELS)),
                ],
            }
    (F7 / "F7_safeguard_A_maps.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
