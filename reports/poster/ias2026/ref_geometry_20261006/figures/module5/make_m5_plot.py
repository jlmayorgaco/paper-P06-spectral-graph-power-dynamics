# M5 plot (design v3): diverging lollipop of the signed eigenvalues of delta D_H, grouped by PLL site.
# source: ../../../one_mode_20261005/generated/data/TABLE_02_FINITE_DELAY_RANK_LAW.csv (lambda_min, lambda_max; 60 rows)
# usage: python make_m5_plot.py W_px H_px   (reference px; saved 1:1 as m5_signed_eigs.pdf)
import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, font_manager as fm
FD = "C:/Users/walla/AppData/Roaming/MiKTeX/fonts/opentype/public/fira/"
for f in ["FiraSans-Regular", "FiraSans-Medium", "FiraSans-SemiBold"]:
    fm.fontManager.addfont(FD + f + ".otf")
W_PX, H_PX = float(sys.argv[1]), float(sys.argv[2])
PXW, PXH = 0.787597/25.4, 0.787299/25.4
rcParams.update({"font.family": "Fira Sans", "font.weight": "medium", "font.size": 19, "pdf.fonttype": 42,
                 "axes.linewidth": 1.4, "xtick.major.width": 0, "ytick.major.width": 1.2,
                 "xtick.major.size": 0, "ytick.major.size": 5, "ytick.direction": "out"})
d = pd.read_csv("../../../one_mode_20261005/generated/data/TABLE_02_FINITE_DELAY_RANK_LAW.csv")
d = d.sort_values(["bus", "design", "frequency_Hz"]).reset_index(drop=True)
TEAL, GOLD, GOLDT, TEXT, GRID, MUTED = "#2C7A70", "#E0A91E", "#9A6F08", "#1C2B27", "#E3EBE7", "#4F5D59"
fig = plt.figure(figsize=(W_PX*PXW, H_PX*PXH))
ax = fig.add_axes([0.155, 0.15, 0.842, 0.84])
sites = sorted(d.bus.unique()); step = 1.0; gap = 0.45
xs = []; centres = []; pos = 0.0
for s in sites:
    n = (d.bus == s).sum()
    for k in range(n): xs.append(pos + k*step)
    centres.append(pos + (n-1)*step/2); pos += (n-1)*step + 1.0 + gap*2
d["x"] = xs
ax.set_yscale("symlog", linthresh=0.1, linscale=0.4)
ax.axhline(0, color=MUTED, lw=1.6, zorder=2)
for _, r in d.iterrows():
    ax.plot([r.x, r.x], [0, r.lambda_max], color=TEAL, lw=3.0, solid_capstyle="round", zorder=3)
    ax.plot([r.x, r.x], [0, r.lambda_min], color=GOLD, lw=3.0, solid_capstyle="round", zorder=3)
ax.scatter(d.x, d.lambda_max, s=46, color=TEAL, edgecolors="white", linewidths=0.8, zorder=4)
ax.scatter(d.x, d.lambda_min, s=46, color=GOLD, edgecolors="white", linewidths=0.8, zorder=4)
ax.set_yticks([-1000, -10, 0, 10])
ax.set_yticklabels(["\u22121000", "\u221210", "0", "10"], color=MUTED)
ax.yaxis.grid(True, color=GRID, linewidth=1.0, zorder=0); ax.set_axisbelow(True)
ax.set_xticks(centres); ax.set_xticklabels([str(s) for s in sites], color=TEXT, fontweight="semibold")
ax.tick_params(axis="x", pad=3)
for sp in ["top", "right", "bottom"]: ax.spines[sp].set_visible(False)
ax.set_ylim(-4500, 450); ax.set_xlim(-0.8, max(xs)+0.8)
ax.text(0.025, 1.0, "60 cases: 10 PLL sites (bus) × 3 frequencies × 2 designs", transform=ax.transAxes, color=MUTED, fontsize=17.5, va="top", ha="left")
y0 = ax.get_position().y0; hh = ax.get_position().height
fig.text(0.022, y0 + hh*0.74, "positive", rotation=90, color=TEAL, fontsize=19, fontweight="semibold", va="center", ha="center")
fig.text(0.022, y0 + hh*0.25, "negative", rotation=90, color=GOLDT, fontsize=19, fontweight="semibold", va="center", ha="center")
fig.savefig("m5_signed_eigs.pdf", transparent=True)
print(d.positive_count.sum(), d.negative_count.sum(), len(d))
