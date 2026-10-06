# M5 plot (design v2, section 5): signed eigenvalues of delta D_H for the 60 frozen cases.
# source: ../../../one_mode_20261005/generated/data/TABLE_02_FINITE_DELAY_RANK_LAW.csv (lambda_min, lambda_max)
# usage: python make_m5_plot.py W_px H_px   (reference px; saved 1:1 as m5_signed_eigs.pdf)
import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, font_manager as fm
FD = "C:/Users/walla/AppData/Roaming/MiKTeX/fonts/opentype/public/fira/"
for f in ["FiraSans-Regular", "FiraSans-Medium", "FiraSans-SemiBold"]:
    fm.fontManager.addfont(FD + f + ".otf")
W_PX, H_PX = float(sys.argv[1]), float(sys.argv[2])
PXW, PXH = 0.787597/25.4, 0.787299/25.4   # inches per reference px
rcParams.update({"font.family": "Fira Sans", "font.weight": "medium", "font.size": 19, "pdf.fonttype": 42,
                 "axes.linewidth": 1.4, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
                 "xtick.major.size": 5, "ytick.major.size": 5, "xtick.direction": "out", "ytick.direction": "out"})
d = pd.read_csv("../../../one_mode_20261005/generated/data/TABLE_02_FINITE_DELAY_RANK_LAW.csv")
d = d.sort_values(["frequency_Hz", "design", "bus"]).reset_index(drop=True)
TEAL, GOLD, TEXT, GRID, MUTED = "#2C7A70", "#D9A520", "#1C2B27", "#E3EBE7", "#4F5D59"
fig = plt.figure(figsize=(W_PX*PXW, H_PX*PXH))
ax = fig.add_axes([0.17, 0.19, 0.825, 0.79])
x = []; pos = 0; centres = []
for f in [0.5, 5.0, 10.0]:
    start = pos
    for design in ["baseline", "analytic"]:
        for k in range(10):
            x.append(pos); pos += 1
        pos += 0.6
    centres.append((start + pos - 1.6)/2)
    pos += 1.0
d["x"] = x
ax.scatter(d.x, d.lambda_max, s=62, color=TEAL, edgecolors="white", linewidths=1.0, zorder=3)
ax.scatter(d.x, d.lambda_min, s=62, color=GOLD, edgecolors="white", linewidths=1.0, zorder=3)
ax.set_yscale("symlog", linthresh=0.1, linscale=0.5)
ax.set_yticks([-1000, -10, 0, 10])
ax.set_yticklabels(["\u22121000", "\u221210", "0", "10"], color=MUTED)
ax.yaxis.grid(True, color=GRID, linewidth=0.9, zorder=0); ax.set_axisbelow(True)
ax.axhline(0, color=MUTED, lw=1.4, zorder=1)
ax.set_xticks(centres); ax.set_xticklabels(["0.5 Hz", "5 Hz", "10 Hz"], color=MUTED)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.set_ylim(-5000, 3000); ax.set_xlim(-1, max(x)+1)
ax.text(0.012, 0.93, "positive direction", transform=ax.transAxes, color=TEAL, fontsize=19, fontweight="semibold", va="top", ha="left")
ax.text(0.012, 0.06, "negative direction", transform=ax.transAxes, color="#A97E10", fontsize=19, fontweight="semibold", va="bottom", ha="left")
ax.set_ylabel("eigenvalues\nof \u03b4$D_H$", fontsize=17.5, color=TEXT, labelpad=3, linespacing=1.0)
fig.savefig("m5_signed_eigs.pdf", transparent=True)
print(d.positive_count.sum(), d.negative_count.sum(), len(d))
