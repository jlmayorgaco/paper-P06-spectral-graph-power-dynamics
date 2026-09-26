"""Gate 3, DECLARED POST-HOC, descriptive: modes behind the 12-plant hyperedges.

The preregistered secondary check found hyperedges of H_RHP outside the four
candidates (kappa_RHP = 6, 9, 11 at g = 0, 0.25, 1 with k = 1). For every one of
those hyperedges this script lists its RHP eigenvalues, their class
(APERIODIC / IN_BAND / BAND_EDGE_ARTIFACT / OSCILLATORY, G1 thresholds) and the
device group carrying most of the participation. It changes no parameter.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _overnight import choose_workers, pin_blas_threads
from G1_f8_rhp_followup import anatomy
from G3_ieee68 import OUT, solve
from ibr_cycles.diagnosis.composability import parse_hypergraph_label


def task(job):
    g, members = job
    case = solve(members, g, 1.0)
    values, right = np.linalg.eig(case.system.A)
    left = np.linalg.inv(right).T
    return [
        {"g": g, "hyperedge": "+".join(map(str, members)), "size": len(members), **m}
        for m in anatomy(values, right, left, case.system.labels, members)
    ]


def main(argv) -> int:
    pin_blas_threads()
    summary = json.loads((OUT / "G3_summary.json").read_text(encoding="utf-8"))[
        "family12"
    ]
    jobs = [
        (float(g), tuple(sorted(e)))
        for g, row in summary.items()
        for e in parse_hypergraph_label(row["label"])
    ]
    with Pool(choose_workers(), initializer=pin_blas_threads) as pool:
        rows = [r for rs in pool.map(task, jobs, chunksize=4) for r in rs]
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "G3c_family_anatomy_posthoc.csv", index=False)
    out = {"label": "DECLARED POST-HOC, descriptive"}
    for g, grp in frame.groupby("g"):
        out[str(g)] = {
            "hyperedges": int(grp.hyperedge.nunique()),
            "modes": int(len(grp)),
            "classes": grp["class"].value_counts().to_dict(),
            "dominant_group": grp.dominant_group.value_counts().to_dict(),
            "f_hz_range": [float(grp.f_hz.min()), float(grp.f_hz.max())],
            "re_range": [float(grp.re.min()), float(grp.re.max())],
            "top_states": grp.top_state.str.rsplit("_", n=1)
            .str[0]
            .value_counts()
            .head(5)
            .to_dict(),
            "plant_9_in_every_edge": bool(
                all("9" in e.split("+") for e in grp.hyperedge.unique())
            ),
        }
    (OUT / "G3c_family_anatomy_posthoc.json").write_text(
        json.dumps(out, indent=1, default=str), encoding="utf-8"
    )
    print(json.dumps(out, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
