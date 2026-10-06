"""Assemble required tables/figures from campaign outputs (no new computation of results)."""
import sys, pathlib, glob, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = CAMP / 'derived'; F = CAMP / 'figures'
# ---------- TABLE09 sparse candidates (engine v2) with gain effort
rows = []
for f in sorted(glob.glob(str(D / 'TABLE09b_final_*.csv'))) + sorted(glob.glob(str(D / 'TABLE09_S0_*.csv'))) + sorted(glob.glob(str(D / 'TABLE09c_scan_*.csv'))):
    t = pd.read_csv(f)
    for r in t.itertuples():
        kp = np.array([float(x) for x in r.Kp.split()]); ki = np.array([float(x) for x in r.Ki.split()])
        eff = float(np.hypot(np.log(kp / KPB), np.log(ki / KIB)).__pow__(2).sum() ** .5)
        changed = int(np.sum((np.abs(np.log(kp / KPB)) > 0.02) | (np.abs(np.log(ki / KIB)) > 0.02)))
        rows.append(dict(label=r.label, sites=r.sites, kind=r.kind, success=r.success, N_unstable=r.N_unstable, N_margin=r.N_margin_excl_gauge, worst_re_catalog=r.worst_re_catalog,
                         gain_effort_l2_logK=eff, sites_changed_gt2pct=changed, iterations=r.iterations))
S = pd.DataFrame(rows).drop_duplicates('label'); S.to_csv(D / 'TABLE09_sparse_candidates.csv', index=False); print(S.round(4).to_string(index=False))
# ---------- TABLE10 root counts: continuation + delay + repaired designs
cnt = pd.read_csv(D / 'TABLE_D01c_root_counts.csv'); cnt.insert(0, 'case', 'baseline gains, uniform delay')
extra = S[['label', 'N_unstable', 'N_margin']].rename(columns={'label': 'case', 'N_margin': 'N_margin_excl_gauge'}); extra['tau_ms'] = 44
pd.concat([cnt, extra], ignore_index=True).to_csv(D / 'TABLE10_root_counts.csv', index=False)
# ---------- FIG08 sparse repair root migration (worst catalogued real part vs iteration)
fig, ax = plt.subplots(figsize=(7.6, 4.6)); cols = {'S0_ki': '#B3363A', 'S0_kp': '#C99A20', 'S0_joint': '#7a3b3b', 'phys4_joint': '#5b6b73', 'core4_joint': '#176494', 'all10_joint': '#003B2D', 'uniform_kp': '#2a8f6b'}
for name, c in cols.items():
    fn = D / f'TABLE09a_S0_history_{name.split("_")[1]}.csv' if name.startswith('S0') else D / f'TABLE09b_hist_{name}.csv'
    if not fn.exists(): continue
    h = pd.read_csv(fn); h = h[h.label == name] if 'label' in h else h
    ax.plot(h.iteration, h.worst_re, '-o', ms=3, color=c, lw=1.3, label=name)
ax.axhline(-SIGMA, color='k', ls='--', lw=1); ax.set_xlabel('trust-region iteration (exact spectrum recomputed each step)'); ax.set_ylabel('worst catalogued Re(lambda) [1/s]'); ax.set_ylim(-0.2, 1.1); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(F / 'FIG08_sparse_repair_root_migration.png', dpi=200); fig.savefig(F / 'FIG08_sparse_repair_root_migration.pdf')
# ---------- TABLE11 + FIG09 forced response
T = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(CAMP / 'raw' / 'forced' / 'forced_*_f*.csv')))]); T = T[~T.label.str.startswith('test')]; T.to_csv(D / 'TABLE11_forced_response.csv', index=False); print(T[['label', 'status', 't_end_s', 'peak_pll_err_deg', 'peak_pll_err_deg_first10s', 'peak_busfreq_dev_Hz_first10s', 't_pll_err_gt_5deg_s']].round(4).to_string(index=False))
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), sharey=False)
for ax, case in zip(axs, ['allSG', 'base40', 'S44']):
    for tag, c in (('f0.7', '#5b6b73'), ('f1.0', '#B3363A'), ('f1.3', '#176494')):
        fn = CAMP / 'raw' / 'forced' / f'trace_{case}_{tag}.csv'
        if fn.exists(): tr = pd.read_csv(fn); ax.semilogy(tr.t, tr.pll_err_deg.clip(lower=1e-6) if case != 'allSG' else tr.busfreq_dev_Hz.clip(lower=1e-7), color=c, lw=1, label={'f0.7': '0.7 f*', 'f1.0': 'f*', 'f1.3': '1.3 f*'}[tag])
    ax.set_title({'allSG': 'all-SG: bus-frequency dev [Hz]', 'base40': '40 ms baseline: PLL error [deg]', 'S44': '44 ms, no retune: PLL error [deg]'}[case], fontsize=9); ax.set_xlabel('time [s]'); ax.grid(alpha=.2, which='both'); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F / 'FIG09_forced_response_allSG_vs_GFL_vs_repaired.png', dpi=200); fig.savefig(F / 'FIG09_forced_response_allSG_vs_GFL_vs_repaired.pdf')
# ---------- TABLE12 events
B = pd.read_csv(ROOT / 'research_gold' / 'checks' / 'T9_FBK_base40.csv').rename(columns={'design': 'label'}); B['label'] = 'EV_baseline_40ms'
E = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(CAMP / 'raw' / 'events' / 'events_*.csv')))] + [B], ignore_index=True); E.to_csv(D / 'TABLE12_nonlinear_events.csv', index=False)
print(E.groupby(E.columns[0]).apply(lambda g: pd.Series(dict(passed=int(g['pass'].astype(str).str.lower().eq('true').sum()), n=len(g), worstF=g.F.max(), worstR=g.R.max(), minslack=g.slack.min()))).to_string())

G = E.groupby('label').apply(lambda g: pd.Series(dict(passed=int(g['pass'].astype(str).str.lower().eq('true').sum()), F=g.F.max(), slack=g.slack.min()))).reset_index()
fig, ax = plt.subplots(figsize=(7.6, 4.2)); y = np.arange(len(G))
ax.barh(y, G.passed, color=['#003B2D' if v == 5 else '#C99A20' for v in G.passed]); ax.set_yticks(y); ax.set_yticklabels(G.label.str.replace('EV_', ''), fontsize=8); ax.set_xlim(0, 5.5)
for yi, (v, f_, sl) in enumerate(zip(G.passed, G.F, G.slack)): ax.text(v + .05, yi, f'{int(v)}/5  max|df|={f_:.3f} Hz  min SG slack={sl:.4f}', va='center', fontsize=7)
ax.set_xlabel('frozen events passed (of 5)'); fig.tight_layout(); fig.savefig(F / 'FIG10_nonlinear_event_summary.png', dpi=200); fig.savefig(F / 'FIG10_nonlinear_event_summary.pdf')
import shutil
for a, b in [('FIG_B01_allSG_spectrum', 'FIG01_allSG_vs_GFL_spectrum'), ('FIG_C01_mode_migration', 'FIG02_replacement_mode_migration'), ('FIG_D01_delay_root_locus', 'FIG03_delay_continuation')]:
    for ext in ('png', 'pdf'): shutil.copy(F / f'{a}.{ext}', F / f'{b}.{ext}')
for a, b in [('TABLE04a_baseline_GFL_modes_40ms', 'TABLE04_GFL_modes'), ('TABLE_D01_delay_continuation', 'TABLE05_delay_continuation'), ('TABLE_F01_feedback_edges', 'TABLE06_feedback_edges'), ('TABLE_F02_pair_cycles', 'TABLE07_pair_cycles'), ('TABLE03_allSG_modes', 'TABLE03_allSG_modes')]:
    if a != b: shutil.copy(D / f'{a}.csv', D / f'{b}.csv')
print(G.to_string())
