"""Phase I: preregistered support S0 = {30,33,36,37}, Ki-only / Kp-only / joint at tau = 44 ms. Everything is recorded, success or not."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from repair import *
S0 = [30, 33, 36, 37]; H = []; F = []
kinds = sys.argv[1:] or ['ki', 'kp', 'joint']
for kind in kinds:
    p, h, f = repair(S0, kind, label=f'S0_{kind}'); H += h; F.append(f)
    pd.DataFrame(H).to_csv(CAMP / 'derived' / f'TABLE09a_S0_history_{"_".join(kinds)}.csv', index=False); pd.DataFrame(F).to_csv(CAMP / 'derived' / f'TABLE09_S0_{"_".join(kinds)}.csv', index=False)
    np.save(CAMP / 'raw' / f'S0_{kind}_p.npy', p)
