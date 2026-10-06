# Hasse lattice of the 16 subsets of H4={30,33,35,37}, coloured by alpha_perp at the frozen policy P4.
import numpy as np, pandas as pd, matplotlib, itertools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import rcParams
rcParams["font.family"] = "Segoe UI"
D = "../../../research/bnd_h4_mechanism/results/20260925T204017_549c3c07_h4_crossmode_v3/"
df = pd.read_csv(D + "derived/TX4_PROPER_SUBSET_CLOSURE.csv")
ALPHA_H4 = 0.12700646782830968  # derived/../claims/tx4_summary.json core.alpha_P4
H4 = [30, 33, 35, 37]
val = {}
for _, r in df.iterrows():
    key = frozenset() if r.subset == "BASE" else frozenset(int(x) for x in str(r.subset).split("+"))
    val[key] = float(r.alpha_perp_at_P4)
    assert bool(r.stable) and r.alpha_perp_at_P4 < 0
val[frozenset(H4)] = ALPHA_H4
assert len(val) == 16
print("proper subsets:", sum(1 for k in val if k != frozenset(H4)), "all stable:", all(v < 0 for k, v in val.items() if k != frozenset(H4)))
print("range stable:", min(v for k, v in val.items() if k != frozenset(H4)), max(v for k, v in val.items() if k != frozenset(H4)))
W, Hh = 234 * 0.787597 / 25.4, 160 * 0.787299 / 25.4
fig = plt.figure(figsize=(W, Hh)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 234); ax.set_ylim(160, 0); ax.axis("off")
levels = {c: [frozenset(s) for s in itertools.combinations(H4, c)] for c in range(5)}
ypos = {4: 15, 3: 49, 2: 83, 1: 117, 0: 145}
nw, nh = 36, 25
cen = {}
for c, subs in levels.items():
    m = len(subs)
    for i, s in enumerate(subs):
        cen[s] = (234 * (i + 0.5) / m, ypos[c])
for s in val:
    for t in val:
        if len(t) == len(s) + 1 and s < t:
            (x0, y0), (x1, y1) = cen[s], cen[t]
            ax.plot([x0, x1], [y0, y1], color="#B9CCD2", lw=1.0, zorder=1)
lo, hi = -0.22, -0.14
for s, (x, y) in cen.items():
    v = val[s]
    if s == frozenset(H4):
        fc, tc, ec = "#C73E3E", "white", "#7a1f1f"
    else:
        u = np.clip((v - lo) / (hi - lo), 0, 1)       # 0 = most damped, 1 = least damped
        base = np.array([0x2C, 0x7A, 0x70]) / 255.0
        fc = tuple(1 - (1 - base) * (0.95 - 0.55 * u)); tc, ec = "#18302A", "#2C7A70"
    ax.add_patch(FancyBboxPatch((x - nw / 2, y - nh / 2), nw, nh, boxstyle="round,pad=0,rounding_size=3", fc=fc, ec=ec, lw=1.0, zorder=3))
    for j, b in enumerate(H4):
        px = x - 12 + j * 8
        ax.scatter(px, y - 7.0, s=22, c=("#18302A" if (b in s and s != frozenset(H4)) else ("white" if b in s else "none")),
                   ec=("white" if s == frozenset(H4) else "#18302A"), lw=0.7, zorder=4)
    txt = ("+%.3f" % v) if v > 0 else ("%.2f" % v).replace("-", "−")
    ax.text(x, y + 5.5, txt, ha="center", va="center", color=tc, fontsize=15.5, fontweight="bold", zorder=5)
fig.savefig("hasse.pdf", transparent=True)
