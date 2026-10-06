"""Central figure of the 20261005 poster: spectrum AND events for four retuning strategies (IEEE-39 design A).
Curves: research_gold/checks/T10_tracked_curves.csv (tracked rightmost catalogued root).
Event markers: EVENTS below, transcribed from research_gold/checks/T5*.csv and T9_*.csv (events meeting every guard / events run)."""
import pathlib, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = pathlib.Path(__file__).resolve().parent
CH = HERE.parents[3] / 'research_gold' / 'checks'
d = pd.read_csv(CH / 'T10_tracked_curves.csv')
BLUE, GREEN, GOLD, RED, GREY = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73'
# (strategy, delay ms): (passed, run)
EVENTS = {('node_4p91', 42): (5, 5), ('node_4p91', 44): (2, 5), ('node_1p04', 44): (5, 5), ('node_1p04', 46): (5, 5), ('node_1p04', 48): (4, 5),
          ('lowfreq', 44): 'aborted', ('none', 44): 'aborted'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 20, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 2.0})
fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.2, 6.3), sharex=True, gridspec_kw=dict(height_ratios=[3.0, 1.55], hspace=0.07))
ax.axhspan(0, 4, color=RED, alpha=.08, lw=0); ax.axhline(0, color=RED, lw=1.6)
sty = {'none': ('no retune', GREY), 'lowfreq': ('low-frequency rule', GOLD),
       'node_4p91': ('protect the 4.91 Hz mode', BLUE), 'node_1p04': ('protect the 1.04 Hz mode', GREEN)}
for k, (lab, c) in sty.items():
    g = d[d.strategy == k]
    if k == 'node_4p91': g = g[g.tau_ms <= 48.99]
    ax.plot(g.tau_ms, g.alpha, color=c, lw=3.8, label=lab, zorder=3)
ax.text(50.6, 0.25, 'unstable', color=RED, fontsize=17, va='bottom', ha='right')
ax.set_xlim(40, 50.8); ax.set_ylim(-0.45, 4.0); ax.set_ylabel(r'rightmost root $\Re\lambda$ [s$^{-1}$]')
ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.16, 1.0), fontsize=16, handlelength=1.3); ax.grid(alpha=.2)
rows = ['node_1p04', 'node_4p91', 'lowfreq', 'none']
for r, k in enumerate(rows):
    bx.text(39.85, r, sty[k][0].replace('protect the ', '').replace(' mode', ''), ha='right', va='center', fontsize=15, color=sty[k][1], fontweight='bold')
    bx.axhline(r, color='#d9e2dd', lw=1, zorder=0)
for (k, t), v in EVENTS.items():
    r = rows.index(k); ok = isinstance(v, tuple) and v[0] == v[1]; txt = f'{v[0]}/{v[1]}' if isinstance(v, tuple) else v
    bx.text(t, r, txt, ha='center', va='center', fontsize=16, fontweight='bold', color='white',
            bbox=dict(boxstyle='round,pad=0.28', fc=GREEN if ok else RED, ec='none'))
bx.set_ylim(-0.6, len(rows) - 0.4); bx.set_yticks([]); bx.spines['left'].set_visible(False)
bx.set_xlabel('PLL delay at all ten sites [ms]')
bx.text(50.75, 3.0, 'events meeting' + chr(10) + 'all limits', ha='right', va='center', fontsize=14, color=GREY)
OUT = HERE / 'generated' / 'figures'
fig.savefig(OUT / 'central.pdf', bbox_inches='tight'); fig.savefig(OUT / 'central.png', bbox_inches='tight', dpi=220); print('ok')
