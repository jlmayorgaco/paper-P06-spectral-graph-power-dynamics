"""Gate 3, DECLARED POST-HOC, descriptive: how far is the 68-bus map from a boundary?

The preregistered G3 map found H_RHP = EMPTY at every node. This script reports
the spectral abscissa of the worst of the 16 portfolios at every node, and which
portfolio attains it. It changes no parameter, and it is not a test: it
distinguishes "far from any boundary" from "barely composable". It uses the same
exact assembly as G3_ieee68.py.
"""

from __future__ import annotations

import json
import math
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _overnight import choose_workers, pin_blas_threads
from G3_ieee68 import G_GRID, K_GRID, OUT, SUBSETS, ZERO, N, parts


def main(argv) -> int:
    pin_blas_threads()
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        table = {m: (a, dg, dk) for m, a, dg, dk in pool.map(parts, SUBSETS)}
    rows = []
    for i in range(N):
        for j in range(N):
            g, k = float(G_GRID[i]), float(K_GRID[j])
            worst, arg, freq = -np.inf, "", float("nan")
            for m in SUBSETS[1:]:
                a0, dg, dk = table[m]
                v = np.linalg.eigvals(a0 + g * dg + (k - 1.0) * dk)
                v = v[np.abs(v) > ZERO]
                c = v[np.argmax(v.real)]
                if c.real > worst:
                    worst, arg, freq = (
                        float(c.real),
                        "+".join(map(str, m)),
                        abs(c.imag) / (2 * math.pi),
                    )
            a0, dg, dk = table[()]
            vb = np.linalg.eigvals(a0 + g * dg + (k - 1.0) * dk)
            rows.append(
                {
                    "g": g,
                    "k": k,
                    "worst_abscissa": worst,
                    "worst_portfolio": arg,
                    "worst_freq_hz": freq,
                    "base_abscissa": float(vb[np.abs(vb) > ZERO].real.max()),
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "G3b_margin_posthoc.csv", index=False)
    summary = {
        "label": "DECLARED POST-HOC, descriptive",
        "worst_abscissa_max": float(frame.worst_abscissa.max()),
        "worst_abscissa_at": frame.loc[
            frame.worst_abscissa.idxmax(),
            ["g", "k", "worst_portfolio", "worst_freq_hz"],
        ].to_dict(),
        "worst_abscissa_min": float(frame.worst_abscissa.min()),
        "worst_portfolio_counts": frame.worst_portfolio.value_counts().to_dict(),
        "worst_freq_range_hz": [
            float(frame.worst_freq_hz.min()),
            float(frame.worst_freq_hz.max()),
        ],
        "base_abscissa_range": [
            float(frame.base_abscissa.min()),
            float(frame.base_abscissa.max()),
        ],
    }
    (OUT / "G3b_margin_posthoc.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
