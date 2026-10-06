"""Phase J: comparators under identical engine, bounds and tau = 44 ms. usage: run_J.py NAME  (NAME in SUPPORTS)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from repair import *
ALL = list(range(30, 40))
SUPPORTS = {
    # name: (sites, kind, tied)
    'uniform_joint': (ALL, 'joint', True), 'uniform_ki': (ALL, 'ki', True), 'uniform_kp': (ALL, 'kp', True),
    'all10_joint': (ALL, 'joint', False), 'all10_ki': (ALL, 'ki', False), 'all10_kp': (ALL, 'kp', False),
    # rules frozen before these runs (see AMENDMENT_02): sensitivity top-4 at the leading root, physical top-2 pairs, core top-2 cycles at the leading root
    'sens4_joint': ([30, 37, 39, 38], 'joint', False), 'sens4_ki': ([30, 37, 39, 38], 'ki', False),
    'phys4_joint': ([33, 34, 35, 36], 'joint', False), 'phys4_ki': ([33, 34, 35, 36], 'ki', False),
    'core4_joint': ([30, 35, 36, 37], 'joint', False), 'core4_ki': ([30, 35, 36, 37], 'ki', False),
}
if __name__ == '__main__':
    for name in sys.argv[1:]:
        sites, kind, tied = SUPPORTS[name]; p, h, f = repair(sites, kind, label=name, tied=tied); f['tied'] = tied
        pd.DataFrame(h).to_csv(CAMP / 'derived' / f'TABLE09b_hist_{name}.csv', index=False); pd.DataFrame([f]).to_csv(CAMP / 'derived' / f'TABLE09b_final_{name}.csv', index=False); np.save(CAMP / 'raw' / f'J_{name}_p.npy', p)
