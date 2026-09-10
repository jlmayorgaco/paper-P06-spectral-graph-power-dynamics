"""Band-limited composability versus full small-signal stability, both benchmarks.

H_Gamma counts right-half-plane modes inside the 0.3-1.5 Hz window only. This
measures how often a point with H = EMPTY (no subset disturbs Gamma) still has a
right-half-plane mode OUTSIDE the window for some subset: F7 frozen maps (from
their stored per-subset RHP counts) and F12 Kundur nodes (recomputed).
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _overnight import choose_workers, pin_blas_threads
from F7_report import points_file
from F12_kundur import G_GRID, MAPS, SUBSETS, Y_GRID, solve, theta_of


def _rhp(task):
    m, i, j = task
    th = theta_of(m, G_GRID[i], Y_GRID[j])
    out = {"map": m, "i": i, "j": j}
    worst = 0
    for s in SUBSETS:
        v = np.linalg.eigvals(solve(s, th).system.A)
        n = int(np.count_nonzero((v.real > 0) & (np.abs(v) > 1e-3) & (v.imag >= 0)))
        out["rhp_" + ("+".join(map(str, s)) or "BASE")] = n
        worst = max(worst, n)
    out["any_rhp"] = worst > 0
    return out


def main() -> int:
    report = {}
    for name in ("F7A", "F7B", "F7C"):
        p = pd.read_csv(points_file(name), low_memory=False)
        ok = p[p.label == "EMPTY"]
        cols = [c for c in p.columns if c.startswith("rhp_") and c != "rhp_BASE"]
        report[name] = {
            "EMPTY_points": int(len(ok)),
            "EMPTY_with_rhp_outside_band": int((ok[cols].max(axis=1) > 0).sum()),
        }
    nodes = pd.read_csv(RESULTS / "F12" / "F12_nodes.csv")
    empty = nodes[nodes.label == "EMPTY"]
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        rows = pool.map(
            _rhp, [(r.map, int(r.i), int(r.j)) for r in empty.itertuples()], chunksize=4
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(RESULTS / "F12" / "F12_empty_rhp.csv", index=False)
    for m in MAPS:
        f = frame[frame["map"] == m]
        report[m] = {
            "EMPTY_points": int(len(f)),
            "EMPTY_with_rhp_outside_band": int(f.any_rhp.sum()),
        }
    (RESULTS / "F12" / "F12_outside_band.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
