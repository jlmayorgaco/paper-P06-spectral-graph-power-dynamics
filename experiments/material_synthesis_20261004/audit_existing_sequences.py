"""Read-only combinatorial audit of frozen TX4 same-policy spectral labels.

No dynamic simulation, gain optimization, or new stability labels are computed.
Paths below mean chains of equilibrium labels, not validated switching dynamics.
"""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
SOURCE = ROOT / 'results' / 'TX4_V9_GLOBAL_VS_EM_TRUTH.csv'
rows = list(csv.DictReader(SOURCE.open(encoding='utf-8-sig', newline='')))
data = {}
for row in rows:
    sites = frozenset() if row['portfolio'] == 'BASE' else frozenset(map(int, row['portfolio'].split('+')))
    assert sites not in data
    assert len(sites) == int(row['cardinality'])
    assert (float(row['alpha_global']) < 0) == (row['global_status'] == 'STABLE')
    data[sites] = row
nodes = sorted(set().union(*data))
assert len(data) == 2 ** len(nodes)
order = sorted(data, key=lambda s: (len(s), tuple(sorted(s))))
paths, one_path, records = {}, {}, []
for s in order:
    stable = data[s]['global_status'] == 'STABLE'
    predecessors = [s - {i} for i in sorted(s)]
    paths[s] = (1 if stable else 0) if not s else (sum(paths[p] for p in predecessors) if stable else 0)
    one_path[s] = []
    if s and paths[s]:
        chosen = next(i for i in sorted(s) if paths[s - {i}])
        one_path[s] = one_path[s - {chosen}] + [chosen]
    records.append({
        'portfolio': data[s]['portfolio'], 'cardinality': len(s),
        'alpha_global': data[s]['alpha_global'], 'global_status': data[s]['global_status'],
        'stable_equilibrium_order_count': paths[s], 'all_order_count': math.factorial(len(s)),
        'reachable_through_stable_equilibria': paths[s] > 0,
        'all_equilibrium_orders_stable': paths[s] == math.factorial(len(s)),
        'one_stable_equilibrium_order': '+'.join(map(str, one_path[s])),
    })
restorations = []
for s in order:
    if data[s]['global_status'] != 'UNSTABLE':
        continue
    for i in nodes:
        t = s | {i}
        if i not in s and data[t]['global_status'] == 'STABLE':
            restorations.append({'unstable_portfolio': data[s]['portfolio'], 'added_GFL_bus': i,
                                 'stable_portfolio': data[t]['portfolio'],
                                 'alpha_before': float(data[s]['alpha_global']),
                                 'alpha_after': float(data[t]['alpha_global'])})
stable_sets = [r for r in records if r['global_status'] == 'STABLE']
conditional = [r for r in stable_sets if r['reachable_through_stable_equilibria'] and not r['all_equilibrium_orders_stable']]
for row in conditional:
    target = frozenset(map(int, row['portfolio'].split('+')))
    direct_count = sum(all(data[frozenset(p[:k])]['global_status'] == 'STABLE'
                           for k in range(len(p) + 1)) for p in itertools.permutations(target))
    assert direct_count == row['stable_equilibrium_order_count']

def chain(sequence):
    return [{'portfolio': data[frozenset(sequence[:k])]['portfolio'],
             'alpha_global': float(data[frozenset(sequence[:k])]['alpha_global']),
             'global_status': data[frozenset(sequence[:k])]['global_status']}
            for k in range(len(sequence) + 1)]

summary = {
    'source': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'policy': '(g,k,t,h)=(0.03625,1.425,1.5,1), frozen TX4 P4',
    'status': 'NUMERICALLY_VALIDATED',
    'scope': 'Exact combinatorics of existing full-spectrum numerical labels; no switching/transient certification.',
    'nodes': nodes, 'portfolio_count': len(data), 'stable_portfolio_count': len(stable_sets),
    'stable_and_reachable_count': sum(r['reachable_through_stable_equilibria'] for r in stable_sets),
    'stable_but_unreachable_count': sum(not r['reachable_through_stable_equilibria'] for r in stable_sets),
    'stable_all_orders_count': sum(r['all_equilibrium_orders_stable'] for r in stable_sets),
    'stable_order_dependent_count': len(conditional),
    'unstable_to_stable_one_addition_count': len(restorations),
    'smallest_absolute_alpha': min(abs(float(r['alpha_global'])) for r in rows),
    'order_dependent_examples': sorted(conditional, key=lambda r: (-r['cardinality'], r['portfolio']))[:3],
    'restoration_examples': restorations[:5],
    'h4_same_final_portfolio_two_orders': {
        'stable_equilibria_order': chain([34, 30, 33, 35, 37]),
        'unstable_intermediate_order': chain([30, 33, 35, 37, 34]),
        'interpretation': 'Stationary spectral labels only; transitions between equilibria are not simulated.'
    },
    'verification': 'Dynamic programming counts agree with independent direct permutation enumeration for both order-dependent targets.',
}
with (OUT / 'EXISTING_SEQUENCE_AUDIT.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
(OUT / 'EXISTING_SEQUENCE_AUDIT.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, indent=2))
