# M8 plot A: residual eps_i(n) vs hop radius n for m = 1..5 target modes, frozen G_c.
# source: research/nhop_ieee39_20261006/derived/E1_residuals.csv (graph=frozen, column eps) and E1_nmin.csv
# usage: python make_m8_plotA.py W_px H_px
from common import *
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
W, H = float(sys.argv[1]), float(sys.argv[2])
d = pd.read_csv(DER + "E1_residuals.csv"); d = d[d.graph == "frozen"]
fig = plt.figure(figsize=(W*PXW, H*PXH))
ax = fig.add_axes([0.135, 0.30, 0.66, 0.66])
FLOOR = 1e-15
col = {1: GREEN, 2: TEAL, 3: NAVY, 4: GOLD, 5: RED}
for m in [1, 2, 3, 4, 5]:
    s = d[d.m == m].sort_values("n")
    y = np.maximum(s.eps.values, FLOOR)
    ax.plot(s.n, y, color=col[m], lw=3.4, zorder=3 if m in (2, 4, 5) else 2, solid_capstyle="round", ls=(0, (1.6, 1.4)) if m == 2 else "-")
    ax.scatter(s.n, y, s=70, color=col[m], edgecolors="white", linewidths=1.0, zorder=4)
ax.axhline(1e-8, color=MUTED, lw=1.4, ls=(0, (4, 3)), zorder=1)
ax.text(2.15, 4e-8, "exact ≤ 10$^{-8}$", color=MUTED, fontsize=17.5, va="bottom", ha="left")
ax.set_yscale("log"); ax.set_ylim(3e-17, 1.2e-1); ax.set_xlim(-0.2, 3.1)
ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["n=0","n=1","n=2","n=3"]); ax.set_yticks([1e-15, 1e-8, 1e-2])
ax.set_yticklabels(["10$^{-15}$", "10$^{-8}$", "10$^{-2}$"], color=MUTED)
ax.tick_params(axis="y", which="minor", length=0); ax.tick_params(axis="x", colors=MUTED, pad=3)
ax.yaxis.grid(True, which="major", color=GRID, lw=0.9, zorder=0); ax.set_axisbelow(True)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.set_ylabel("residual ε", color=TEXT, fontsize=17.5, labelpad=2)
kw = dict(fontsize=17.5, fontweight="semibold", ha="center", va="bottom")
ax.text(0.02, 2e-13, "m=1", color=GREEN, **kw)
ax.text(1, 2e-13, "m=2, 3", color=NAVY, **kw)
ax.text(2, 2e-13, "m=4", color="#A97E10", **kw)
ax.text(3.2, 4.8e-4, "m=5:\nnot reached", color=RED, fontsize=17.5, fontweight="semibold", ha="left", va="center", linespacing=0.95)
fig.savefig("m8_plotA.pdf", transparent=True)
for m in range(1, 6): print(m, d[d.m == m].sort_values("n").eps.values)
