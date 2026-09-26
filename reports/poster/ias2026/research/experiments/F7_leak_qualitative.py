"""F7 check - do the qualitative F7 phenomena survive the regulator leak?

F7_leak_sensitivity.py shows that boundary LOCATIONS move with the leak. This
re-traces the four signature pure-policy lines and the coarse F7A hypergraph
count at w = 0.02, 0.05 and 0.20 on the fast path. Run from the research root.
"""

import sys

sys.path.insert(0, "src")
sys.path.insert(0, "experiments")  # noqa: E702
import numpy as np
import pandas as pd

import _f7_common as C
from ibr_cycles.diagnosis.composability import hypergraph_order, parse_hypergraph_label

d = pd.read_csv("results/F7/F7_leak_sensitivity.csv")
for w in ("0.02", "0.2"):
    bad = d[d.label != d[f"label_w{w}"]]
    print(
        "w",
        w,
        "disagreements by map/g:",
        bad.groupby("map")
        .g.describe()[["count", "min", "50%", "max"]]
        .round(3)
        .to_dict("index"),
    )
gs = np.round(np.concatenate([np.arange(0, 0.1, 0.0025), np.arange(0.1, 0.6, 0.01)]), 4)
lines = {
    "F7B t=0.852": dict(k=1.0, t=0.852, h=1.0),
    "F7A k=1.30": dict(k=1.30, t=1.5, h=1.0),
    "F7C h=0.5": dict(k=1.0, t=1.0, h=0.5),
    "F7B t=0.68": dict(k=1.0, t=0.68, h=1.0),
}
for w in (0.02, 0.05, 0.2):
    C.LEAK = w
    C._FAST = None
    for name, fx in lines.items():
        labs = [C.evaluate(C.Theta(g=g, **fx))["label"] for g in gs]
        order = [hypergraph_order(parse_hypergraph_label(lab)) for lab in labs]
        kap = ["inf" if k == -1 else str(k) for k in order]
        seq = [k for i, k in enumerate(kap) if i == 0 or k != kap[i - 1]]
        hs = [x for i, x in enumerate(labs) if i == 0 or x != labs[i - 1]]
        print(f"w={w} {name}: kappa {' > '.join(seq)}  | distinct H {len(set(labs))}")
    # distinct H over a coarse F7A-like grid
    lab = set()
    for g in (0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5):
        for k in np.arange(0.5, 2.25, 0.05):
            lab.add(C.evaluate(C.Theta(g=g, k=k, t=1.5))["label"])
    print(f"w={w} distinct H on coarse F7A grid: {len(lab - {'BASE_UNSTABLE'})}")
