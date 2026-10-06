"""T10: tracked rightmost catalogued root vs delay for four strategies (design A). Continuation only; counts come from T3/T6/T7/T8."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
from t1_spillover_ieee39 import transport, roots
fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real); mid = min(roots, key=lambda z: abs(z.imag - 6.544))
rows = []
for name, rule, lam in [('none', 'none', None), ('lowfreq', 'lowfreq', None), ('node_4p91', 'exact', fast), ('node_1p04', 'exact', mid)]:
    cur = list(roots)
    for h in np.arange(0, 0.01001, 0.00025):
        p = P0.copy()
        for i in range(10):
            if rule == 'exact': p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], lam, h)
            elif rule == 'lowfreq': p[10 + i] = P0[10 + i] + h * P0[20 + i]
        m = Model(p, TAU + h); new = []
        for z in cur:
            try: new.append(m.refine(z)[0])
            except Exception: new.append(z)
        cur = new; w = max(cur, key=lambda z: z.real)
        rows.append(dict(strategy=name, tau_ms=40 + h * 1e3, alpha=w.real, f_Hz=w.imag / 2 / np.pi, Kp=p[10], Ki=p[20]))
pd.DataFrame(rows).to_csv('T10_tracked_curves.csv', index=False); print('ok')
