# M9 plot (design v2, section 5): delay root-locus, real part of the PLL-family roots vs uniform delay, full IEEE-39.
# source: ../../../../../../research/feedback_cycle_fullmodel_20261005/derived/TABLE_D01_delay_continuation.csv (branch, tau_ms, real, PLL_part)
#         ../../../../../../research/feedback_cycle_fullmodel_20261005/derived/TABLE_D01b_delay_crossings.csv (tau_cross_margin_ms)
# PLL family = branches whose mean PLL_part > 0.9 (R0-R5, R21, R38, R55); margin = -0.05 1/s (FINAL_FEEDBACK_CYCLE_REPORT.md section 11)
# usage: python make_m9_plot.py W_px H_px   (reference px; saved 1:1 as m9_delay_locus.pdf)
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
                 "axes.linewidth": 1.4, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
                 "xtick.major.size": 5, "ytick.major.size": 5, "xtick.direction": "out", "ytick.direction": "out"})
R = "../../../../../../research/feedback_cycle_fullmodel_20261005/derived/"
d = pd.read_csv(R + "TABLE_D01_delay_continuation.csv")
cr = pd.read_csv(R + "TABLE_D01b_delay_crossings.csv")
GREEN, NAVY, RED, GRID, MUTED, TEXT, GOLD = "#03534A", "#1E4F8A", "#B83A3A", "#E3EBE7", "#4F5D59", "#1C2B27", "#D9A520"
pll = d.groupby("branch").PLL_part.mean()
pll = pll[pll > 0.9].index
p = d[d.branch.isin(pll)]
tau_c = float(cr[cr.branch == "R0"].tau_cross_margin_ms.iloc[0])
fig = plt.figure(figsize=(W_PX*PXW, H_PX*PXH))
ax = fig.add_axes([0.14, 0.27, 0.83, 0.71])
ax.axhspan(-0.05, 3.3, color="#F6DEDE", zorder=0, lw=0)
for b, g in p.groupby("branch"):
    ax.plot(g.tau_ms, g.real, color=NAVY, lw=1.8, alpha=0.6, zorder=2, solid_capstyle="round")
env = p.groupby("tau_ms").real.max()
ax.plot(env.index, env.values, color=GREEN, lw=4.0, zorder=3, solid_capstyle="round")
ax.axhline(-0.05, color=RED, lw=2.0, zorder=2)
ax.axvline(tau_c, color=GOLD, lw=2.4, ls=(0, (3, 2)), zorder=1)
ax.plot([tau_c], [-0.05], "o", ms=15, mfc=GOLD, mec=GREEN, mew=2.0, zorder=5)
ax.set_xlim(40, 50); ax.set_ylim(-3.1, 3.3)
ax.set_xticks([40, 42, 44, 46, 48, 50]); ax.set_yticks([-2, 0, 2])
ax.set_yticklabels(["\u22122", "0", "2"])
ax.tick_params(colors=MUTED, labelsize=17.5)
ax.yaxis.grid(True, color=GRID, linewidth=0.9, zorder=0); ax.set_axisbelow(True)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.text(40.3, 3.0, "unstable", color=RED, fontsize=17.5, fontweight="semibold", ha="left", va="top", zorder=6)
ax.text(40.2, 0.35, "stable", color=GREEN, fontsize=18, fontweight="semibold", ha="left", va="bottom") if False else None
ax.set_xlabel("uniform delay \u03c4 (ms)", fontsize=17.5, color=TEXT, labelpad=2)
ax.set_ylabel("Re \u03bb (1/s)", fontsize=17.5, color=TEXT, labelpad=2)
fig.savefig("m9_delay_locus.pdf", transparent=True)
print(tau_c, len(pll), env.loc[44.0])
