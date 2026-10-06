# M8 plot B: roots left beyond the -0.05 margin after an exact assignment, frozen G_c.
# source: research/nhop_ieee39_20261006/derived/E5_spillover.csv (graph=frozen; N_margin; variants absolute, delta); baseline 10 roots beyond margin
#         (REPORT_E1_E2_E5.md: "Roots beyond -0.05 left (baseline 10, unstable 8)").
# usage: python make_m8_plotB.py W_px H_px
from common import *
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
W, H = float(sys.argv[1]), float(sys.argv[2])
d = pd.read_csv(DER + "E5_spillover.csv"); d = d[d.graph == "frozen"]
a = d[d.variant == "absolute"].sort_values(["m", "n"]).reset_index(drop=True)
dl = d[d.variant == "delta"].sort_values(["m", "n"]).reset_index(drop=True)
assert (a.m.values == dl.m.values).all() and (a.n.values == dl.n.values).all() and len(a) == 12
x = []; pos = 0; prev = None
for m in a.m:
    if prev is not None and m != prev: pos += 0.7
    x.append(pos); pos += 1; prev = m
x = np.array(x)
fig = plt.figure(figsize=(W*PXW, H*PXH))
ax = fig.add_axes([0.01, 0.25, 0.98, 0.72])
YM = 22
cols = [TEAL if v == 0 else RED for v in a.N_margin]
ax.bar(x, np.where(a.N_margin == 0, 0.5, a.N_margin), width=0.72, color=cols, zorder=3)
ax.scatter(x, dl.N_margin, s=120, marker="D", color=NAVY, edgecolors="white", linewidths=1.4, zorder=5)
for xi, v in zip(x, a.N_margin):
    ax.text(xi, v + 0.5 if v else 1.1, str(int(v)) if v else "none", color=RED if v else TEAL, fontsize=17.5, fontweight="bold", ha="center", va="bottom")
ax.axhline(10, color=MUTED, lw=1.6, ls=(0, (4, 3)), zorder=2)
ax.text(x[-1] + 0.62, 9.4, "no retune: 10", color=MUTED, fontsize=17.5, fontweight="semibold", ha="right", va="top")
ax.set_ylim(0, YM); ax.set_xlim(-0.6, x[-1] + 0.7)
ax.set_yticks([]); ax.set_xticks(x); ax.set_xticklabels([str(int(n)) for n in a.n], color=MUTED)
ax.tick_params(axis="x", length=0, pad=2)
for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
for m in range(1, 5):
    xs = x[a.m.values == m]
    ax.annotate(f"$m{{=}}{m}$", xy=(xs.mean(), 0), xycoords=("data", "axes fraction"), xytext=(0, -25), textcoords="offset points",
                color=TEXT, fontsize=17.5, fontweight="semibold", ha="center", va="top", annotation_clip=False)
ax.text(-0.5, 21.5, "x-axis: n = hop radius", color=MUTED, fontsize=17.5, ha="left", va="top") if False else None
fig.savefig("m8_plotB.pdf", transparent=True)
print(list(zip(a.m, a.n, a.N_margin, dl.N_margin)), (a.N_margin > 0).sum())
