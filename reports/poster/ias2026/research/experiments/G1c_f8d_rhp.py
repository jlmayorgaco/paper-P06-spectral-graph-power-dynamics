"""Journal gate 1, completion: the F8D AVR-speed path recomputed under Gamma_RHP.

G1_whole_rhp.py recomputed the F8 factorials. The F8D continuation (every
surviving AVR slowed, efd' -> beta efd', at the six F8 points, D = 0 and 2) is
recomputed here with the same configurations, counting eigenvalues in the whole
right half plane (|s| > 1e-3). The frozen band labels are read from
results/F8/F8D_avr_blend.csv and are not modified.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import SUBSETS, Theta
from _overnight import pin_blas_threads
from F8_service_attribution import solve_config
from G1_whole_rhp import label_of, rhp_modes

OUT = RESULTS / "G1"
WORKERS = 12


def job(task):
    pname, theta_d, beta, damping = task
    cfg = {
        "class": "A",
        "name": "blend",
        "services": {
            "A1_inertia": 1.0,
            "A2_damping": damping,
            "A3_flux": 1.0,
            "A4_avr": beta,
            "A5_pss": 1.0,
        },
    }
    theta = Theta(**theta_d)
    counts = {}
    for members in SUBSETS:
        v = np.linalg.eigvals(solve_config(members, theta, cfg).system.A)
        if not members and v[np.abs(v) > 1e-3].real.max() >= 0:
            return {
                "point": pname,
                "beta": beta,
                "D": damping,
                "H_RHP": "BASE_UNSTABLE",
            }
        counts[frozenset(members)] = int(rhp_modes(v).size)
    return {"point": pname, "beta": beta, "D": damping, "H_RHP": label_of(counts)}


def main(argv) -> int:
    pin_blas_threads()
    pts = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    frozen = pd.read_csv(RESULTS / "F8" / "F8D_avr_blend.csv")
    tasks = [
        (
            r.point,
            {k: pts[r.point][k] for k in ("g", "k", "t", "h")},
            float(r.beta),
            r.D,
        )
        for r in frozen.itertuples()
    ]
    with Pool(WORKERS, initializer=pin_blas_threads) as pool:
        rows = pd.DataFrame(pool.map(job, tasks, chunksize=2))
    rows = rows.merge(
        frozen[["point", "beta", "D", "label"]], on=["point", "beta", "D"]
    )
    rows.to_csv(OUT / "G1_f8d_rhp.csv", index=False)
    seqs = {}
    for (p, d), g in rows.groupby(["point", "D"]):
        g = g.sort_values("beta", ascending=False)
        empty = g[g.H_RHP == "EMPTY"].beta
        seqs[f"{p}|D={d:g}"] = {
            "H_RHP_empty_at_beta_max": float(empty.max()) if len(empty) else None,
            "band_empty_at_beta_max": float(g[g.label == "EMPTY"].beta.max())
            if (g.label == "EMPTY").any()
            else None,
            "differences": int((g.H_RHP != g.label).sum()),
        }
    summary = {
        "cases": int(len(rows)),
        "H_RHP_equals_band": float((rows.H_RHP == rows.label).mean()),
        "per_point": seqs,
    }
    (OUT / "G1_f8d_rhp.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
