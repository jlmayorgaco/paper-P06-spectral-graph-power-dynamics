"""POST-HOC sensitivity (not preregistered): H with Gamma = the whole right half plane.

The preregistered F12 test uses the 0.3-1.5 Hz window. On Kundur every pair of
replacements carries an aperiodic (zero-frequency) right-half-plane mode that the
window excludes by construction, so this recomputes the hypergraph with every
RHP mode counted, on the Kundur nodes (direct path) and on the frozen F7 maps
(stored per-subset RHP counts), and asks whether policy still moves it.
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import LABELS as F7_LABELS
from _f7_common import SUBSETS as F7_SUBSETS
from _overnight import choose_workers, pin_blas_threads
from F7_report import points_file
from F12_kundur import G_GRID, MAPS, SUBSETS, Y_GRID, solve, theta_of
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    incompatibility_hypergraph,
)


def _node(task):
    m, i, j = task
    th = theta_of(m, G_GRID[i], Y_GRID[j])
    counts = {}
    for s in SUBSETS:
        v = np.linalg.eigvals(solve(s, th).system.A)
        counts[frozenset(s)] = int(
            np.count_nonzero((v.real > 0) & (np.abs(v) > 1e-3) & (v.imag >= 0))
        )
    if counts[frozenset()] > 0:
        return {"map": m, "i": i, "j": j, "label_full": "BASE_UNSTABLE"}
    return {
        "map": m,
        "i": i,
        "j": j,
        "label_full": hypergraph_label(incompatibility_hypergraph(counts)),
    }


def main() -> int:
    report = {}
    nodes = pd.read_csv(RESULTS / "F12" / "F12_nodes.csv")
    ok = nodes[~nodes.label.isin(["BASE_UNSTABLE", "INFEASIBLE"])]
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        rows = pool.map(
            _node, [(r.map, int(r.i), int(r.j)) for r in ok.itertuples()], chunksize=8
        )
    full = ok[["map", "i", "j", "g", "label"]].merge(
        pd.DataFrame(rows), on=["map", "i", "j"]
    )
    full.to_csv(RESULTS / "F12" / "F12_full_rhp.csv", index=False)
    for m, y in MAPS.items():
        f = full[full["map"] == m]
        f = f.merge(nodes[nodes["map"] == m][["i", "j", y]], on=["i", "j"])
        lines = f.groupby("j").label_full.nunique()
        report[m] = {
            "nodes": int(len(f)),
            "distinct_H_full": int(f.label_full.nunique()),
            "H_full_counts": f.label_full.value_counts().to_dict(),
            "policy_lines_with_2plus_H_full": int((lines > 1).sum()),
            "policy_lines": int(len(lines)),
            "EMPTY_full": int((f.label_full == "EMPTY").sum()),
            "same_as_band_label": float((f.label_full == f.label).mean()),
        }
    for name in ("F7A", "F7B", "F7C"):
        p = pd.read_csv(points_file(name), low_memory=False)
        p = p[~p.label.isin(["BASE_UNSTABLE", "INFEASIBLE"])]
        labels = []
        for r in p[[f"rhp_{lab}" for lab in F7_LABELS]].to_numpy():
            labels.append(
                hypergraph_label(
                    incompatibility_hypergraph(
                        {
                            frozenset(s): int(n)
                            for s, n in zip(F7_SUBSETS, r, strict=True)
                        }
                    )
                )
            )
        p = p.assign(label_full=labels)
        report[name] = {
            "points": int(len(p)),
            "same_as_band_label": float((p.label_full == p.label).mean()),
            "distinct_H_full": int(p.label_full.nunique()),
            "distinct_H_band": int(p.label.nunique()),
        }
    (RESULTS / "F12" / "F12_full_rhp_sensitivity.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
