"""F8 Q4 - is the surviving fleet's AVR necessary, or one coordinate among others?

Manual excitation (field held) makes the BASE grid aperiodically unstable, so it
is not a viable configuration (F8C). Here the AVR of every surviving machine is
slowed continuously, efd' -> beta efd', from the native AVR (beta = 1) towards
manual excitation, at the F8 points, and H_Gamma is recorded together with the
base-stability limit. Direct path.
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import SUBSETS, Theta, in_gamma
from _overnight import choose_workers, pin_blas_threads
from F8_service_attribution import solve_config
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    incompatibility_hypergraph,
)

OUT = RESULTS / "F8"
BETAS = np.round(np.geomspace(1.0, 0.002, 40), 5)


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
        if not members:
            dyn = v[np.abs(v) > 1e-3]
            base = float(dyn.real.max())
            if base >= 0:
                return {
                    "point": pname,
                    "beta": beta,
                    "D": damping,
                    "label": "BASE_UNSTABLE",
                    "base_abscissa": base,
                }
        counts[frozenset(members)] = int(np.count_nonzero(in_gamma(v)))
    return {
        "point": pname,
        "beta": beta,
        "D": damping,
        "base_abscissa": base,
        "label": hypergraph_label(incompatibility_hypergraph(counts)),
    }


def main() -> int:
    pts = json.loads((OUT / "F8_points.json").read_text(encoding="utf-8"))
    tasks = [
        (n, {k: p[k] for k in ("g", "k", "t", "h")}, float(b), d)
        for n, p in pts.items()
        for b in BETAS
        for d in (0.0, 2.0)
    ]
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        rows = pool.map(job, tasks, chunksize=2)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "F8D_avr_blend.csv", index=False)
    for (p, d), g in frame.groupby(["point", "D"]):
        g = g.sort_values("beta", ascending=False)
        labs = g.label.tolist()
        seq = [
            f"{lab}@{b:g}"
            for i, (lab, b) in enumerate(zip(labs, g.beta, strict=True))
            if i == 0 or lab != labs[i - 1]
        ]
        print(p, "D", d, "|", " > ".join(seq))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
