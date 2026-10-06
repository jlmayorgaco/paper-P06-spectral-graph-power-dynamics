"""T16: catalogued roots in the s-plane at 40 ms (base) and 44 ms for four strategies (continuation, design A)."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
from t1_spillover_ieee39 import transport, roots
fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real); mid = min(roots, key=lambda z: abs(z.imag - 6.544))
rows = [dict(strategy='base40', re=z.real, im=z.imag) for z in roots]
for name, rule, lam in [('none', 'none', None), ('lowfreq', 'lowfreq', None), ('node491', 'exact', fast), ('node104', 'exact', mid)]:
    cur = list(roots)
    for h in np.linspace(0, 0.004, 17)[1:]:
        p = P0.copy()
        for i in range(10):
            if rule == 'exact': p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], lam, h)
            elif rule == 'lowfreq': p[10 + i] = P0[10 + i] + h * P0[20 + i]
        m = Model(p, TAU + h); new = []
        for z in cur:
            try: new.append(m.refine(z)[0])
            except Exception: new.append(z)
        cur = new
    rows += [dict(strategy=name, re=z.real, im=z.imag, origin_re=z0.real, origin_im=z0.imag) for z, z0 in zip(cur, roots)]
pd.DataFrame(rows).to_csv('T16_splane_roots.csv', index=False); print('ok', fast, mid)
