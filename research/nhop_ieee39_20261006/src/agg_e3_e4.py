"""Aggregate raw/runs/*.json -> derived/E3_repair.csv, derived/E4_replacement.csv, figures."""
import json, pathlib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
CAMP = pathlib.Path(__file__).resolve().parents[1]
rows = [json.loads(f.read_text())['final'] for f in sorted((CAMP / 'raw' / 'runs').glob('*.json'))]
df = pd.DataFrame(rows)
df['fd_max_rel_err'] = df.fd_rel_err.apply(lambda d: max(d.values()))
df = df.drop(columns='fd_rel_err')
cols = ['graph', 'n', 'tau_ms', 'rho', 'success', 'N_unstable', 'N_margin', 'iterations', 'worst_catalog_re', 'worst_catalog_f_hz', 'effort', 'n_offdiag_vars', 'fd_max_rel_err', 'max_noninteger_winding', 'runtime_s', 'stop_note', 'label']
df = df[cols].sort_values(['graph', 'tau_ms', 'rho', 'n'])
E3 = df[df.rho == 0.875]; E3.to_csv(CAMP / 'derived' / 'E3_repair.csv', index=False)
E4 = df[df.tau_ms == 44]; E4.to_csv(CAMP / 'derived' / 'E4_replacement.csv', index=False)
# largest verified rho per (graph, n)
LR = []
for (g, n), d in E4.groupby(['graph', 'n']):
    ok = d[d.success]; LR.append(dict(graph=g, n=n, largest_verified_rho=float(ok.rho.max()) if len(ok) else np.nan,
                                      note='largest verified point in the declared sweep, not a certified optimum'))
LR = pd.DataFrame(LR); LR.to_csv(CAMP / 'derived' / 'E4_largest_rho.csv', index=False)
print(df.drop(columns=['label', 'stop_note']).to_string()); print(LR.to_string())
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 18})
W, H = 136 / 25.4, 84 / 25.4
cols_ = {44: 'C0', 48: 'C1', 52: 'C2'}
fig, ax = plt.subplots(figsize=(W, H))
for g, ls, mf in (('frozen', '-', True), ('bridge', '--', False)):
    for tau in (44, 48, 52):
        d = E3[(E3.graph == g) & (E3.tau_ms == tau)].sort_values('n'); off = 0.06 * (1 if g == 'bridge' else -1)
        ax.plot(d.n + off, d.worst_catalog_re, ls, color=cols_[tau], marker='o', mfc=cols_[tau] if mf else 'none', ms=7, label=f'{g} {tau} ms')
        f = d[~d.success]; ax.plot(f.n + off, f.worst_catalog_re, 'x', color='k', ms=13, mew=2)
ax.axhline(-0.05, color='gray', ls=':'); ax.set_yscale('symlog', linthresh=0.1); ax.set_xlabel('hop radius n'); ax.set_ylabel('worst root Re (1/s)')
ax.set_xticks(range(5)); ax.legend(fontsize=8, ncol=2, loc='upper left', title='x = full-spectrum failure', title_fontsize=8)
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'E3_repair.pdf'); plt.close(fig)
fig, ax = plt.subplots(figsize=(W, H))
for g, mk, off in (('frozen', 'o', -0.04), ('bridge', 's', 0.04)):
    d = LR[LR.graph == g]; ax.plot(d.n + off, d.largest_verified_rho, mk + '-', label=g)
ax.set_ylim(0.86, 0.96); ax.set_xlabel('hop radius n'); ax.set_ylabel('largest verified rho'); ax.set_xticks(range(5)); ax.legend(fontsize=14, loc='lower right')
ax.text(0.5, 0.97, 'largest verified point in the declared sweep', transform=ax.transAxes, ha='center', va='top', fontsize=9)
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'E4_replacement.pdf'); plt.close(fig)
