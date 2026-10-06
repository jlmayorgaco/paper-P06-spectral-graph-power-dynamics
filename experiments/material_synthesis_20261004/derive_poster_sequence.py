"""Exact finite-state widest-path calculation from existing TX4 labels.

No simulation or new spectral classification. Initial equilibrium is included
in each path's minimum decay margin. These are equilibrium sequences, not
certified switching trajectories.
"""
from pathlib import Path
import csv
import hashlib
import itertools
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / 'results' / 'TX4_V9_GLOBAL_VS_EM_TRUTH.csv'
rows = list(csv.DictReader(SOURCE.open(encoding='utf-8-sig', newline='')))
alpha = {frozenset() if r['portfolio'] == 'BASE' else frozenset(map(int, r['portfolio'].split('+'))):
         float(r['alpha_global']) for r in rows}
assert len(alpha) == 512
best, witness = {}, {}
for s in sorted(alpha, key=lambda x: (len(x), tuple(sorted(x)))):
    if not s:
        best[s], witness[s] = -alpha[s], ()
    else:
        candidates = [(min(-alpha[s], best[s - {i}]), witness[s - {i}] + (i,)) for i in s]
        best[s], witness[s] = min(candidates, key=lambda item: (-item[0], item[1]))

target = frozenset([30, 33, 34, 35, 37])
h4 = frozenset([30, 33, 35, 37])
j34 = frozenset([34])
bad = (30, 33, 35, 37, 34)
def path_alphas(order):
    return [alpha[frozenset(order[:k])] for k in range(len(order) + 1)]

enumeration = [(p, min(-a for a in path_alphas(p))) for p in itertools.permutations(sorted(target))]
assert abs(max(b for _, b in enumeration) - best[target]) < 1e-14
assert min(-a for a in path_alphas(witness[target])) == best[target]
summary = {
    'source': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'model': 'TX4 frozen P4, complete transverse spectrum of equilibria after removing rotational/drift modes; not the later delayed model',
    'status': 'NUMERICALLY_VALIDATED',
    'theorem_status': 'EXACT_IDENTITY for the finite optimization defined on this table',
    'objective': 'maximize the minimum of -alpha_global over all prefix equilibria, including BASE',
    'target': sorted(target), 'optimal_order': witness[target],
    'optimal_minimum_decay_margin_s-1': best[target],
    'optimal_order_alpha_s-1': path_alphas(witness[target]),
    'bad_order': bad, 'bad_order_alpha_s-1': path_alphas(bad),
    'bad_minimum_decay_margin_s-1': min(-a for a in path_alphas(bad)),
    'orders_total': len(enumeration),
    'orders_all_prefixes_stable': sum(b > 0 for _, b in enumeration),
    'orders_attaining_optimum_within_1e-12': sum(abs(b - best[target]) < 1e-12 for _, b in enumeration),
    'contextual_bus34_effect': {
        'alpha_base': alpha[frozenset()], 'alpha_34_only': alpha[j34],
        'alpha_H4': alpha[h4], 'alpha_H4_plus_34': alpha[target],
        'delta_alpha_at_base': alpha[j34] - alpha[frozenset()],
        'delta_alpha_at_H4': alpha[target] - alpha[h4],
        'two_block_interaction': alpha[target] - alpha[h4] - alpha[j34] + alpha[frozenset()],
        'scope': 'Finite marginal of spectral abscissa; mode identity may switch. H4 is one action block here, not an irreducible fifth-order interaction.'
    },
    'optimality_upper_bound': 'Every path includes BASE, so its minimum margin cannot exceed -alpha(BASE); the witness attains this upper bound.',
    'verification': 'Dynamic program independently agrees with all 120 permutations.',
    'limitations': ['fixed gains and binary replacement candidates',
                    'no numerical eigenvalue enclosures or continuous replacement interpolation',
                    'no nonlinear transition, event, voltage, current, or RoCoF validation',
                    'not an optimum for MW, gain tuning, or replacement fraction',
                    'TX4 SG has D=0 and no governor; these are transverse synchronization results, not frequency-security certificates'],
}
(HERE / 'POSTER_SEQUENCE_DESIGN.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(10.6, 5.6), constrained_layout=True)
ax.axhspan(0, .16, color='#fff0ee', zorder=0)
ax.axhline(0, color='#6b7280', linewidth=1)
ax.plot(range(6), path_alphas(bad), 'o-', color='#c64c3c', lw=2.3,
        label='Bus 34 last: unstable intermediate equilibrium')
ax.plot(range(6), path_alphas(witness[target]), 's-', color='#176b80', lw=2.3,
        label='Best minimum-margin order: 30 → 33 → 34 → 35 → 37')
ax.annotate('H4: +0.1270 s$^{-1}$', xy=(4, alpha[frozenset([30,33,35,37])]),
            xytext=(2.1, .095), arrowprops={'arrowstyle':'->','color':'#c64c3c'}, color='#9e3025')
ax.annotate('Same final design\n−0.2050 s$^{-1}$', xy=(5, alpha[target]),
            xytext=(3.6, -.105), arrowprops={'arrowstyle':'->','color':'#374151'}, color='#374151')
ax.set(title='More GFL can restore transverse stability — replacement order matters',
       xlabel='Number of completed SG → GFL substitutions',
       ylabel='Complete transverse spectral abscissa α (s$^{-1}$)', ylim=(-.25,.16), xlim=(-.12,5.15))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.grid(axis='y', alpha=.18)
ax.spines[['right','top']].set_visible(False)
ax.legend(loc='lower left', fontsize=9)
fig.suptitle('IEEE-39 · fixed TX4 P4 policy · equilibrium spectra, with rotation/drift removed', fontsize=10, color='#555b65')
fig.savefig(HERE / 'FIG_POSTER_SEQUENCE_CONCEPT.png', dpi=220)
plt.close(fig)
print(json.dumps(summary, indent=2))
