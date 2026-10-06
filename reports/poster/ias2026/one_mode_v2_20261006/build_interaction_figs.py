"""Poster figures for the interaction counterexample and repair, from archived
experiments/interaction_decision_20261004 tables (no recomputation)."""
import json, pathlib
import pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve().parents[4] / 'experiments' / 'interaction_decision_20261004'
OUT = pathlib.Path(__file__).resolve().parent / 'generated' / 'figures'
enc = pd.read_csv(R / 'TABLE_15_ROOT_ENCLOSURES.csv').set_index('design')
it = pd.read_csv(R / 'TABLE_06_COMPENSATION_ITERATION.csv')
BLUE, GREEN, GOLD, RED, GREY = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 19, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 1.6})
x = lambda k: enc.loc[k, 'root_real']
# --- A: finite decision reversal
fig, ax = plt.subplots(figsize=(11.4, 4.6))
ax.axvspan(-0.05, -0.036, color=RED, alpha=.10, lw=0)
ax.axvline(-0.05, color=RED, ls='--', lw=2.4)
ax.text(-0.0492, 3.62, 'margin violated', color=RED, fontsize=19, va='center')
rows = [('anchor', 'Anchor design', GREY), ('single30', 'Retune bus 30 only', BLUE),
        ('single37', 'Retune bus 37 only', BLUE), ('joint', 'Both retunes together', RED),
        ('corrected', 'Interaction-corrected', GREEN)]
for k, (key, lab, c) in enumerate(rows):
    ax.scatter(x(key), k, s=330, color=c, zorder=4)
ax.annotate('', xy=(x('joint'), 3 - .0), xytext=(x('anchor'), 3), arrowprops=dict(arrowstyle='->', lw=2.6, color=RED))
ax.text(-0.0600, 2.62, r'$\mathcal{I}=+0.0243\ \mathrm{s}^{-1}$', ha='center', color=RED, fontsize=21, fontweight='bold')
ax.set_ylim(4.6, -.6); ax.set_yticks(range(5)); ax.set_yticklabels([r[1] for r in rows], fontweight='bold')
for t,(_,_,c) in zip(ax.get_yticklabels(),rows): t.set_color(c)
ax.set_xlim(-0.0675, -0.036)
ax.set_xlabel(r'Real part of the enclosed PLL pole $\Re\lambda$  [s$^{-1}$]')
ax.grid(axis='x', alpha=.2)
fig.savefig(OUT / 'decision.png', bbox_inches='tight', dpi=250); fig.savefig(OUT / 'decision.pdf', bbox_inches='tight'); plt.close(fig)
# --- B: fixed-point repair
fig, ax = plt.subplots(figsize=(10.2, 4.2))
ax.axhline(-0.05, color=RED, ls='--', lw=2.4); ax.axhline(-0.06, color=GREEN, ls=':', lw=2.2)
ax.plot(it.k, it.joint_real, 'o-', color=BLUE, lw=3, ms=11)
ax.text(10.1, -0.0485, 'required margin', color=RED, ha='right', va='bottom', fontsize=18)
ax.text(10.1, -0.0612, 'target $-0.060$', color=GREEN, ha='right', va='top', fontsize=18)
ax.set_xlabel('Fixed-point update $k$'); ax.set_ylabel(r'Joint $\Re\lambda$  [s$^{-1}$]')
ax.set_xticks(range(0, 11, 2)); ax.grid(alpha=.2)
fig.savefig(OUT / 'repair.png', bbox_inches='tight', dpi=250); fig.savefig(OUT / 'repair.pdf', bbox_inches='tight'); plt.close(fig)
print('ok')
