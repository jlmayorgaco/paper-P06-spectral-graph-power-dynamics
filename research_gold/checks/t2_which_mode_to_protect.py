"""T2: finite all-site delay transport on IEEE-39 (40 ms -> 40+h ms at every PLL).
Strategies: none (fixed gains) / protect a slow mode / protect the least-damped ~4.9 Hz mode.
Tracks the catalogued physical roots by continuation (NOT a full-spectrum count)."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS","1"); os.environ.setdefault("MKL_NUM_THREADS","1")
import numpy as np, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU, ROOTS
from t1_spillover_ieee39 import transport, roots
print('Kp0', P0[10], 'Ki0', P0[20])
slow = min(roots, key=lambda z: z.real if z.imag < 1 else 9)       # placeholder, replaced below
slow = [z for z in roots if z.imag < 1][0]; fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real)
rows = []
for name, lam in [('none', None), ('protect_slow', slow), ('protect_fast', fast)]:
    cur = list(roots)
    for h in np.arange(0, 0.01001, 0.00025):
        p = P0.copy(); tau = TAU + h
        if lam is not None:
            for i in range(10): p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], lam, h)
        m = Model(p, tau); new = []
        for z in cur:
            try: new.append(m.refine(z)[0])
            except Exception: new.append(complex('nan'))
        cur = [n if n == n else z for n, z in zip(new, cur)]
        worst = max(cur, key=lambda z: z.real)
        rows.append(dict(strategy=name, h_ms=h * 1e3, alpha=worst.real, f_Hz=worst.imag / 2 / np.pi, Kp=p[10], Ki=p[20],
                         prot_drift=(min(abs(z - lam) for z in cur) if lam is not None else np.nan)))
df = pd.DataFrame(rows); df.to_csv('T2_which_mode_to_protect.csv', index=False)
for n, g in df.groupby('strategy'):
    print(n); print(g[g.h_ms.isin([0, 2, 4, 6, 8, 10])].round(5).to_string(index=False))
for lam, n in [(slow, 'slow'), (fast, 'fast')]:
    a, w = lam.real, lam.imag
    print(n, lam, 'Ki-sign budget h* [ms]:', 1e3 * (np.arctan2(w * P0[20], abs(lam) ** 2 * P0[10] + a * P0[20])) / w)
