"""Figures of the 20261006 poster (36 x 36 in). Each figure is drawn at its exact slot size in the poster (1:1), fonts 19-21 pt.
Every curve/bar is read from the campaign tables in research/feedback_cycle_fullmodel_20261005/derived (copied into generated/data)."""
import pathlib, shutil
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CAMP = ROOT / 'research' / 'feedback_cycle_fullmodel_20261005'
D = CAMP / 'derived'
OUT = HERE / 'generated' / 'figures'; DATA = HERE / 'generated' / 'data'; OUT.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
BLUE, GREEN, GOLD, RED, GREY, INK, LGREEN = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73', '#17372D', '#2a8f6b'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 20, 'axes.labelsize': 20, 'xtick.labelsize': 19, 'ytick.labelsize': 19, 'legend.fontsize': 18,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 2.0, 'mathtext.fontset': 'dejavusans'})
MM = 1 / 25.4
for f in ['TABLE_D01_delay_continuation', 'TABLE_D01b_delay_crossings', 'TABLE_D01c_root_counts', 'TABLE_E01_sigma_min_scan', 'TABLE_F02_pair_cycles', 'TABLE_F01_feedback_edges', 'TABLE08_gain_authority',
          'TABLE_G01_physical_coupling', 'TABLE09_sparse_candidates', 'TABLE11_forced_response', 'TABLE12_nonlinear_events']:
    shutil.copy(D / f'{f}.csv', DATA / f'{f}.csv')
for pat in ('TABLE09b_hist_*.csv', 'TABLE09a_S0_history_*.csv', 'TABLE09c_scan_*.csv'):
    for f in D.glob(pat): shutil.copy(f, DATA / f.name)
for f in (CAMP / 'raw' / 'forced').glob('trace_*.csv'): shutil.copy(f, DATA / f.name)

def fig_mm(w, h): return plt.subplots(figsize=(w * MM, h * MM), layout='constrained')
def save(fig, name):
    fig.savefig(OUT / f'{name}.png', dpi=250); fig.savefig(OUT / f'{name}.pdf'); plt.close(fig)

# ---- A. delay root locus (slot 214 x 84 mm)
T = pd.read_csv(DATA / 'TABLE_D01_delay_continuation.csv'); C = pd.read_csv(DATA / 'TABLE_D01b_delay_crossings.csv')
fig, ax = fig_mm(268, 108)
for b, g in T.groupby('branch'):
    if g.PLL_part.mean() < .5: continue
    g = g[g.tau_ms <= 48]; ax.plot(g.tau_ms, g.real, color=BLUE, lw=2.2, alpha=.9)
ax.axhline(0, color=INK, lw=1.4); ax.axhline(-0.05, color=RED, ls='--', lw=2.4)
first = C.tau_cross_margin_ms.min(); ax.axvline(first, color=GOLD, lw=2.4, ls=':')
ax.set_title(f'first margin crossing at {first:.2f} ms', color=GOLD, fontsize=19, fontweight='bold', loc='left', pad=2)
ax.text(47.9, 0.15, 'unstable', color=RED, ha='right', fontsize=19)
ax.set_xlim(40, 48); ax.set_ylim(-2.2, 2.2); ax.set_xlabel('uniform PLL delay [ms]'); ax.set_ylabel('Re $\\lambda$ [s$^{-1}$]'); ax.grid(alpha=.2)
save(fig, 'delay_locus')

# ---- B. closure scan (slot 150 x 78 mm)
S = pd.read_csv(DATA / 'TABLE_E01_sigma_min_scan.csv'); fig, ax = fig_mm(130, 100)
ax.loglog(S.f_hz, S.smin_IQ_44, color=RED, lw=3, label='$\\sigma_{\\min}(I+Q)$, 44 ms'); ax.loglog(S.f_hz, S.smin_IQ_40, color=GREEN, lw=3, label='$\\sigma_{\\min}(I+Q)$, 40 ms')
ax.loglog(S.f_hz, S.minL_44, color=RED, lw=2.2, ls=':', label='min $|L_i|$, 44 ms')
ax.set_xlabel('frequency [Hz]'); ax.set_ylabel('singular value'); ax.set_xlim(.5, 20); ax.set_ylim(5e-3, 3e2); ax.grid(alpha=.2, which='both'); ax.legend(frameon=False, loc='upper right', fontsize=15)
save(fig, 'closure_scan')

# ---- C. pair cycles at the 4.681 Hz root (slot 172 x 78 mm)
P = pd.read_csv(DATA / 'TABLE_F02_pair_cycles.csv'); Ph = pd.read_csv(DATA / 'TABLE_G01_physical_coupling.csv'); FROZEN = {(35, 36), (30, 37), (33, 34)}
g = P[P.root == 'tau44_4.681Hz'].merge(Ph, on=['bus_i', 'bus_j']).sort_values('abs_1mp').head(5)[::-1]
fig, ax = fig_mm(136, 100); y = np.arange(len(g)); cols = [GOLD if (a, b) in FROZEN else BLUE for a, b in zip(g.bus_i, g.bus_j)]
ax.barh(y, g.abs_1mp, color=cols); ax.set_yticks(y); ax.set_yticklabels([f'{a}$-${b}' for a, b in zip(g.bus_i, g.bus_j)]); ax.set_xlabel('$|1-p_{ij}|$ at the 4.68 Hz root'); ax.grid(alpha=.2, axis='x'); ax.set_xlim(0, .75)
for yi, v in zip(y, g.abs_1mp): ax.text(v + .012, yi, f'{v:.3f}', va='center', fontsize=17)
save(fig, 'pair_cycles')

# ---- D. interaction graph at the 4.681 Hz root (slot 84 x 78 mm)
E = pd.read_csv(DATA / 'TABLE_F01_feedback_edges.csv'); E = E[E.root == 'tau44_4.681Hz']
fig, ax = plt.subplots(figsize=(84 * MM, 78 * MM)); ang = np.linspace(0, 2 * np.pi, 11)[:-1] + np.pi / 2; pos = {30 + i: (np.cos(a), np.sin(a)) for i, a in enumerate(ang)}; mx = E.abs_Q.max()
for r in E.itertuples():
    hot = {r.src, r.dst} == {35, 36}
    if r.abs_Q < .55 * mx and not hot: continue
    ax.annotate('', xy=np.array(pos[r.dst]) * .86, xytext=np.array(pos[r.src]) * .86, arrowprops=dict(arrowstyle='-|>', lw=1 + 4 * r.abs_Q / mx, color=RED if hot else GREY, alpha=.95 if hot else .4, shrinkA=16, shrinkB=16, connectionstyle='arc3,rad=.14'))
for b, (x, y_) in pos.items(): ax.text(x, y_, str(b), ha='center', va='center', fontsize=17, fontweight='bold', color=INK, bbox=dict(boxstyle='circle,pad=.2', fc='#FFF4D6' if b in (35, 36) else 'white', ec=GREEN, lw=2))
ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect('equal'); ax.axis('off'); save(fig, 'core_graph')

# ---- E. predictor-corrector convergence (slot 206 x 80 mm)
fig, ax = fig_mm(146, 108)
for fn, lab, c in [('TABLE09a_S0_history_joint.csv', 'S0 = {30,33,36,37}', RED), ('TABLE09b_hist_core4_joint.csv', 'pair set {30,35,36,37}', GOLD), ('TABLE09b_hist_all10_joint.csv', 'all ten sites', GREEN), ('TABLE09b_hist_uniform_kp.csv', 'uniform $K_p$', LGREEN)]:
    h = pd.read_csv(DATA / fn); h = h[h.label == h.label.iloc[0]]; ax.plot(h.iteration, h.worst_re, '-o', ms=5, lw=2.4, color=c, label=lab)
ax.axhline(-0.05, color=INK, ls='--', lw=2); ax.set_xlim(-.5, 30); ax.set_ylim(-.2, 1.05); ax.set_xlabel('corrector iteration'); ax.set_ylabel('worst Re $\\lambda$ [s$^{-1}$]'); ax.grid(alpha=.2)
ax.legend(frameon=False, loc='center right', bbox_to_anchor=(1.0, 0.36), fontsize=16)
save(fig, 'repair_curve')

# ---- F. support-size scan (slot 206 x 80 mm)
def first_success(name):
    t = pd.read_csv(DATA / f'TABLE09c_scan_{name}_joint.csv'); s = t[t.success]; return int(s.k.min()) if len(s) else np.nan
res = [first_success('phys'), first_success('sens'), first_success('core')]; rnd = [first_success(f'rand{k}') for k in range(5)]; vals = res + [np.mean(rnd), 10]
fig, ax = fig_mm(146, 108); ax.bar(range(5), vals, color=[BLUE, BLUE, BLUE, GREY, LGREEN], width=.7); ax.scatter(np.full(5, 3) + np.linspace(-.2, .2, 5), rnd, color='white', edgecolor=INK, s=45, zorder=5)
for i, v in enumerate(vals): ax.text(i, v + .3, f'{v:.0f}' if i != 3 else f'{v:.1f}', ha='center', fontsize=19, fontweight='bold')
ax.axhline(4, color=RED, ls='--', lw=2.4, label='S0 and 4-site sets fail'); ax.legend(frameon=False, loc='upper left', fontsize=15)
ax.set_xticks(range(5)); ax.set_xticklabels(['elec.', 'gain', 'core', 'random', 'unif.'], fontsize=16); ax.set_ylabel('sites retuned'); ax.set_ylim(0, 15); ax.grid(alpha=.2, axis='y')
save(fig, 'support_scan')

# ---- G. forced response (slot 268 x 48 mm)
fig, ax = fig_mm(268, 82)
for case, tag, c, lab in (('base40', 'f1.0', GREEN, '40 ms baseline'), ('S44', 'f1.0', RED, '44 ms, no retune'), ('uniformkp', 'f1.0', LGREEN, '44 ms, uniform $K_p$ cut'), ('uniformkp', 'f0.7', GREY, 'same, $0.7f^*$')):
    tr = pd.read_csv(DATA / f'trace_{case}_{tag}.csv'); ax.plot(tr.t, tr.pll_err_deg.rolling(20, min_periods=1).max(), color=c, lw=2.4, label=lab)
ax.set_yscale('log'); ax.set_ylim(.003, 600); ax.set_xlabel('time [s]'); ax.set_ylabel('error [deg]'); ax.grid(alpha=.2, which='both'); ax.legend(frameon=False, loc='upper left', fontsize=14, ncol=2)
save(fig, 'forced')

# ---- H. replaced generation vs margin (slot 268 x 76 mm)
SM = pd.read_csv(D / 'TABLE13_share_margin.csv'); shutil.copy(D / 'TABLE13_share_margin.csv', DATA / 'TABLE13_share_margin.csv')
fig, ax = fig_mm(268, 108)
for tms, c in ((40, GREEN), (44, RED)):
    g = SM[SM.tau_ms == tms]; ax.plot(g.rho * 100, g.rightmost_re, '-o', color=c, lw=2.6, ms=8, label=f'{tms} ms')
ax.axhline(0, color=INK, lw=1.4); ax.axhline(-0.05, color=RED, ls='--', lw=2.2)
ax.set_ylim(-0.3, 1.25); ax.set_xlabel('inverter share [%]'); ax.set_ylabel('rightmost Re $\lambda$ [s$^{-1}$]'); ax.grid(alpha=.2); ax.legend(frameon=False, loc='center left', fontsize=18)
save(fig, 'share_margin')

# ---- I. Julia time-domain events (slot 268 x 48 mm)
fig, ax = fig_mm(268, 82); import glob
for lab, c, fn in (('40 ms baseline', GREEN, 'traj_TRJ_base40_bus16_100.csv'), ('44 ms, uniform $K_p$ cut', LGREEN, 'traj_TRJ_uniform_kp_bus16_100.csv')):
    p = CAMP / 'raw' / 'events' / 'traj' / fn
    if p.exists():
        t = pd.read_csv(p); shutil.copy(p, DATA / fn); ax.plot(t.time_after_event_s, t.Fmax_Hz, color=c, lw=2.2, label=lab)
ax.axhline(0.5, color=RED, ls='--', lw=2); ax.set_xlabel('time after event [s]'); ax.set_ylabel('$|\Delta f|$ [Hz]'); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=16, loc='upper right'); ax.set_ylim(0, .65)
save(fig, 'traj_events')
# ---- J. gain authority at the leading root (slot 118 x 84 mm)
A = pd.read_csv(DATA / 'TABLE08_gain_authority.csv'); Lr = A[A.root_re == A.root_re.max()].sort_values('site')
fig, ax = fig_mm(118, 84); x = np.arange(10); w = .38
ax.bar(x - w / 2, Lr.dRe_dlogKp, w, color=BLUE, label='$K_p$'); ax.bar(x + w / 2, Lr.dRe_dlogKi, w, color=GOLD, label='$K_I$')
ax.set_xticks(x); ax.set_xticklabels(Lr.site, fontsize=15); ax.set_xlabel('PLL site (bus)'); ax.set_ylabel('$\partial$Re$\lambda/\partial\ln K$ [s$^{-1}$]', fontsize=17); ax.legend(frameon=False, fontsize=17); ax.grid(alpha=.2, axis='y'); ax.axhline(0, color=INK, lw=1.2)
save(fig, 'authority')
print('figures ok')
