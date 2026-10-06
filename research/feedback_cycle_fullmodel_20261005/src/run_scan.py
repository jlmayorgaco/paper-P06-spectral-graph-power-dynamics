"""Support-size scan along a frozen site order (AMENDMENT_03); stops at the first successful size. usage: run_scan.py ORDER_NAME KIND"""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from repair import *
order = json.load(open(CAMP / 'derived' / 'ORDERS.json'))[sys.argv[1]]; kind = sys.argv[2] if len(sys.argv) > 2 else 'joint'; rows = []
for k in range(1, 11):
    p, h, f = repair(order[:k], kind, label=f'scan_{sys.argv[1]}_{kind}_k{k}', verbose=False); f['order'] = sys.argv[1]; f['k'] = k; rows.append(f); print(k, order[:k], f['success'], f['worst_re_catalog'], f['N_margin_excl_gauge'], flush=True)
    pd.DataFrame(rows).to_csv(CAMP / 'derived' / f'TABLE09c_scan_{sys.argv[1]}_{kind}.csv', index=False)
    if f['success']: np.save(CAMP / 'raw' / f'scan_{sys.argv[1]}_{kind}_p.npy', p); break
