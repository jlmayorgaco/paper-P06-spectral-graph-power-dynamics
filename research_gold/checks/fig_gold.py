"""Central figure for the recommendation: IEEE-39, all ten PLLs, delay 40 -> 50 ms. Data: T2 (tracked rightmost root) and T3 (full-spectrum counts)."""
import pathlib, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = pathlib.Path(__file__).resolve().parent
t2 = pd.read_csv(H / 'T2_which_mode_to_protect.csv'); t3 = pd.read_csv(H / 'T3_full_spectrum_count.csv')
BLUE, GREEN, GOLD, RED, GREY = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73'
plt.rcParams.update({'font.size': 15, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 1.6})
fig, ax = plt.subplots(figsize=(10.5, 5.6))
ax.axhspan(0, 4, color=RED, alpha=.08, lw=0); ax.axhline(0, color=RED, lw=1.5); ax.axhline(-0.05, color=RED, ls='--', lw=1.3)
lab = {'none': ('no retune', GREY), 'protect_slow': ('hide delay from the 0.05 Hz mode\n(first-order / moment compensation)', GOLD),
       'protect_fast': ('hide delay from the 4.91 Hz mode\n(closed-form transport)', GREEN)}
for k, (name, c) in lab.items():
    g = t2[t2.strategy == k]
    if k == 'protect_fast': g = g[g.h_ms <= 8.99]
    ax.plot(40 + g.h_ms, g.alpha, lw=3.4, color=c, label=name)
    for _, r in t3[t3.strategy == k].iterrows():
        y = np.interp(r.tau_ms, 40 + g.h_ms, g.alpha)
        ax.annotate(f"{int(r.unstable_roots)}", (r.tau_ms, y), textcoords='offset points', xytext=(0, 9), ha='center', color=c, fontsize=13, fontweight='bold')
ax.axvline(48.99, color=GREEN, ls=':', lw=2.2); ax.text(49.05, 2.6, 'phase budget' + chr(10) + 'h* = 8.99 ms', color=GREEN, fontsize=13, va='top')
ax.text(40.1, 3.75, 'unstable', color=RED, fontsize=13, va='top'); ax.text(40.1, -0.33, 'numbers: unstable roots (full-spectrum count)', color=GREY, fontsize=12)
ax.set_xlim(40, 50.9); ax.set_ylim(-0.45, 4); ax.set_xlabel('PLL delay at all ten sites [ms]'); ax.set_ylabel('rightmost tracked root, Re λ [s$^{-1}$]')
ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.10, 0.93), fontsize=13); ax.grid(alpha=.2)
fig.savefig(H.parent / 'FIG_GOLD_one_pll_one_mode.png', dpi=200, bbox_inches='tight'); print('ok')
