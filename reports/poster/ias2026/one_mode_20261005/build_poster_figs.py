"""Poster figures for version 20261004. Values come from archived tables; the only computed curve is the
closed-form floor 2*asin(P/2)/P. Figure sizes are small on purpose so the text stays large on the poster."""
import pathlib, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
R = ROOT / 'experiments' / 'interaction_decision_20261004'
OUT = HERE / 'generated' / 'figures'
BLUE, GREEN, GOLD, RED, GREY, INK = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73', '#17372D'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 21, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 2.0, 'mathtext.fontset': 'dejavusans'})


def save(fig, name):
    fig.savefig(OUT / f'{name}.png', bbox_inches='tight', dpi=250)
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------ s-plane (hero)
enc = pd.read_csv(R / 'TABLE_15_ROOT_ENCLOSURES.csv').set_index('design')
pt = lambda k: (enc.loc[k, 'root_real'], enc.loc[k, 'root_imag'])
fig, ax = plt.subplots(figsize=(9.6, 5.0))
ax.axvspan(-0.05, -0.034, color=RED, alpha=.10, lw=0)
ax.axvline(-0.05, color=RED, ls='--', lw=2.8)
a, s30, s37, j, c, k = (pt(x) for x in ['anchor', 'single30', 'single37', 'joint', 'complex_pair', 'corrected'])


def arrow(p, q, col, ls='-'):
    ax.annotate('', xy=q, xytext=p, arrowprops=dict(arrowstyle='-|>', lw=3.0, color=col, ls=ls, shrinkA=14, shrinkB=14))


arrow(a, s37, BLUE, '--'); arrow(a, s30, BLUE, '--'); arrow(a, j, RED); arrow(j, k, GREEN)
ax.scatter(*c, s=700, facecolor='none', edgecolor=GREEN, lw=3.4, zorder=4)
ax.scatter(*a, s=300, color=GREY, zorder=5)
ax.scatter(*s30, s=260, color=BLUE, zorder=5); ax.scatter(*s37, s=260, color=BLUE, zorder=5)
ax.scatter(*j, s=380, color=RED, zorder=5); ax.scatter(*k, s=340, color=GREEN, zorder=5)
T = lambda x, y, t, col, ha='left', fs=21: ax.text(x, y, t, color=col, fontsize=fs, fontweight='bold', ha=ha, va='center')
T(a[0] + 0.0022, 30.488, 'anchor', GREY)
T(s37[0] + 0.0022, 30.585, 'bus 37 alone', BLUE)
T(s30[0] + 0.0022, 30.885, 'bus 30 alone', BLUE)
T(j[0], 31.075, 'both retunes,' + chr(10) + 'decay-matched', RED, ha='center')
T(k[0], 31.075, 'repaired', GREEN, ha='center')
T(-0.0355, 30.49, 'decay below' + chr(10) + 'requirement', RED, ha='right', fs=19)
ax.set_xticks([-0.065, -0.055, -0.045, -0.035])
ax.set_xlim(-0.0685, -0.0345); ax.set_ylim(30.43, 31.13)
ax.set_xlabel(r'Decay of the 4.9 Hz mode, $\Re\lambda$ [s$^{-1}$]')
ax.set_ylabel(r'Frequency $\Im\lambda$ [rad/s]')
ax.grid(alpha=.2)
save(fig, 'splane')

# ------------------------------------------------------------------ network (30 and 37 marked)
coordinates = (HERE / 'reference_template' / 'network_ieee39.tikz').read_text(encoding='utf-8-sig')
pos = {int(b): (float(y), float(x)) for b, x, y in re.findall(r'\\coordinate \(b(\d+)\) at \(([\d.]+)mm,([\d.]+)mm\)', coordinates)}
branches = pd.read_csv(ROOT / 'reports' / 'experiment_D' / 'inputs' / 'branch.csv')
fig, ax = plt.subplots(figsize=(7.4, 4.6))
for _, e in branches.iterrows():
    p, q = pos[int(e.src_bus)], pos[int(e.dst_bus)]
    ax.plot([p[0], q[0]], [p[1], q[1]], color='#a6b9ad', lw=2.2, zorder=1)
for bus, (x, y) in pos.items():
    hi = bus in (30, 37)
    col = GOLD if hi else (BLUE if bus >= 30 else '#e4ede7')
    ax.scatter(x, y, s=420 if hi else (270 if bus >= 30 else 190), c=col, edgecolor=INK, lw=2.2 if hi else 0.9, zorder=3)
    ax.text(x, y, str(bus), ha='center', va='center', fontsize=11 if not hi else 13, color=INK if (hi or bus < 30) else 'white', fontweight='bold' if hi else 'normal', zorder=4)
ax.set_aspect('equal'); ax.axis('off')
ax.legend(handles=[Line2D([], [], marker='o', color='none', markerfacecolor=GOLD, markeredgecolor=INK, markersize=15, label='buses 30, 37'),
                   Line2D([], [], marker='o', color='none', markerfacecolor=BLUE, markersize=13, label='other inverter sites')],
          loc='upper center', bbox_to_anchor=(.5, 0.02), ncol=2, frameon=False, fontsize=17, handletextpad=.3, columnspacing=1.2)
fig.tight_layout(); save(fig, 'network2')

# ------------------------------------------------------------------ fixed-point repair
it = pd.read_csv(R / 'TABLE_06_COMPENSATION_ITERATION.csv')
fig, ax = plt.subplots(figsize=(7.6, 3.9))
ax.axhspan(-0.05, -0.040, color=RED, alpha=.10, lw=0)
ax.axhline(-0.05, color=RED, ls='--', lw=2.8)
ax.axhline(-0.06, color=GREEN, ls=':', lw=2.6)
ax.plot(np.r_[it.k, 11], np.r_[it.joint_real, -0.06], 'o-', color=BLUE, lw=3.2, ms=9)
ax.text(10.6, -0.0486, 'required', color=RED, ha='right', va='bottom', fontsize=19)
ax.text(10.6, -0.0576, 'target', color=GREEN, ha='right', va='bottom', fontsize=19)
ax.set_xlabel('update $k$'); ax.set_ylabel(r'joint $\Re\lambda$ [s$^{-1}$]')
ax.set_xticks(range(0, 12, 2)); ax.set_ylim(-0.0655, -0.0395); ax.grid(alpha=.2)
save(fig, 'repair')

# ------------------------------------------------------------------ rank-two bars
RK = pd.read_csv(HERE / 'generated' / 'data' / 'TABLE_02_FINITE_DELAY_RANK_LAW.csv')
r = RK[(RK.design == 'analytic') & (RK.frequency_Hz == 5.)]
fig, ax = plt.subplots(figsize=(7.4, 3.9))
ax.axhline(0, c=INK, lw=1.2)
ax.bar(r.bus - .2, r.lambda_max, width=.38, color=BLUE, label='positive')
ax.bar(r.bus + .2, r.lambda_min, width=.38, color=RED, label='negative')
ax.set_yscale('symlog', linthresh=1)
ax.set_xticks(r.bus); ax.set_xlabel('PLL site given +1 ms')
ax.set_ylabel(r'eigenvalue [MW]')
ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(.5, 1.22), fontsize=19)
fig.tight_layout(); save(fig, 'rank_two2')

# ------------------------------------------------------------------ inertia floor (closed form)
P = np.linspace(1e-3, 1.93, 400)
thr = 2 * np.arcsin(P / 2) / P
fig, ax = plt.subplots(figsize=(7.6, 4.5))
ax.fill_between(P, 0.9, thr, color=RED, alpha=.11, lw=0)
ax.plot(P, thr, color=GREEN, lw=4, label='exact floor')
ax.axhline(1.0, color=GOLD, ls='--', lw=3, label='linearised floor')
ax.scatter([1.0], [1.025], s=260, color=RED, zorder=5)
ax.annotate('$m=1.025$:' + chr(10) + 'bus 1 must overshoot', xy=(1.0, 1.025), xytext=(0.08, 1.2), fontsize=20, color=RED,
            arrowprops=dict(arrowstyle='->', lw=2.4, color=RED))
ax.text(1.12, 0.93, 'forbidden region', fontsize=18, color=RED, ha='center', va='bottom')
ax.set_xlim(0, 1.95); ax.set_ylim(0.9, 1.4)
ax.set_xlabel('localised step $P$'); ax.set_ylabel('total inertia $m$')
ax.legend(frameon=False, fontsize=18, loc='upper left', bbox_to_anchor=(0.0, 1.0)); ax.grid(alpha=.2)
save(fig, 'floor2')

# ---- share x delay feasible-region map (panel 7); data: research_gold/checks
CH = ROOT / 'research_gold' / 'checks'
T11 = pd.read_csv(CH / 'T11_share_delay_map.csv')
def ev(name):
    f = CH / f'{name}.csv'
    if not f.exists(): return None
    d = pd.read_csv(f); return f'{int(d["pass"].astype(str).str.lower().eq("true").sum())}/{len(d)}'
EVS = {(0.875, 44): '5/5', (0.875, 48): '4/5', (0.85, 44): ev('T9_B_rule_node_44ms'),
       (0.90, 44): ev('T9_M90_rule_44'), (0.90, 48): ev('T9_M90_rule_48'), (0.80, 44): ev('T9_M80_rule_44')}
shares = [0.80, 0.85, 0.875, 0.90, 0.925, 0.95]; delays = [40, 44, 48]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
for j, rho in enumerate(shares):
    r = T11[np.isclose(T11.rho, rho)].iloc[0]
    for i, d in enumerate(delays):
        if r.base_40ms > 0:
            col, txt = '#B0B0B0', ''
        elif d == 40:
            txt = ev('T9_M90_base_40') if np.isclose(rho, 0.90) else ''
            txt = txt or ''
            col = GOLD if (txt and txt.split('/')[0] != txt.split('/')[1]) else GREEN
        else:
            col, txt = (GREEN if r[f'rule_{d}'] == 0 else RED), EVS.get((rho, d)) or ''
            if col == GREEN and txt and txt.split('/')[0] != txt.split('/')[1]: col = GOLD
        ax.add_patch(plt.Rectangle((j - .46, i - .46), .92, .92, color=col, alpha=.85, lw=0))
        if txt: ax.text(j, i, txt, color='white', fontsize=19, fontweight='bold', ha='center', va='center')
        if d > 40 and r.base_40ms == 0:
            ax.plot([j - .4, j + .4], [i - .36, i - .36], color=RED, lw=5, solid_capstyle='butt')
ax.set_xticks(range(len(shares))); ax.set_xticklabels(['80', '85', '87.5', '90', '92.5', '95'], fontsize=19)
ax.set_yticks(range(3)); ax.set_yticklabels(['40', '44', '48'], fontsize=19)
ax.set_xlim(-.5, 5.5); ax.set_ylim(-.5, 2.5)
ax.set_xlabel('inverter share (%)', fontsize=19); ax.set_ylabel('delay (ms)', fontsize=19)
ax.legend(handles=[Line2D([0], [0], color=GREEN, lw=8, label='rule: 0 roots'), Line2D([0], [0], color=GOLD, lw=8, label='0 roots, events fail'), Line2D([0], [0], color=RED, lw=5, label='no retune: roots'),
                   Line2D([0], [0], color='#B0B0B0', lw=8, label='fails at 40 ms')], frameon=False, fontsize=16, loc='upper center', bbox_to_anchor=(0.5, 1.32), ncol=2, handlelength=1.2, columnspacing=1.0)
save(fig, 'share_delay_map')
print('ok')
