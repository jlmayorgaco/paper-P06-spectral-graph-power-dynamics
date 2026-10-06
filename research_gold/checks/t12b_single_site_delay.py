"""T12b: fairness check for the cooperative K_p-only idea in its own scenario (delay grows at ONE PLL only), IEEE-39 design A."""
import os, sys, pathlib, itertools
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
import model as MM
from model import Model, P0, TAU
from t1_spillover_ieee39 import roots, transport
from t3_full_spectrum_count import count
KPMIN, KPMAX = 7.853981633974483, 125.66370614359172
base = Model(P0, TAU); fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real); kp0 = P0[10:20].copy(); rows = []
for site, hms in [(0, 8), (7, 8), (0, 16), (0, 30)]:
    tau = TAU.copy(); tau[site] += hms * 1e-3
    none = count(Model(P0, tau), -0.05)[0] - 1
    row = dict(site=30 + site, extra_ms=hms, none_beyond_guard=none)
    if none > 0:
        p = P0.copy(); p[10 + site], p[20 + site] = transport(P0[10 + site], P0[20 + site], fast, hms * 1e-3)
        row['local_transport_beyond_guard'] = count(Model(p, tau), -0.05)[0] - 1; row['local_Ki'] = p[20 + site]
        def detK(kp):
            rb = la.solve(fast * base.eye - base.A0, MM.M['Bp'] * kp + MM.M['Bi'] * P0[20:30], check_finite=False)
            return np.linalg.det(np.eye(10) - (base.C @ rb) * np.exp(-fast * tau)[None, :])
        sc = abs(detK(kp0)); best = None
        for i, j in itertools.combinations([q for q in range(10) if q != site], 2):
            def f(x, y):
                k = kp0.copy(); k[i] += x; k[j] += y; return detK(k) / sc
            c0 = f(0, 0); c1 = f(1, 0) - c0; c2 = f(0, 1) - c0; c12 = f(1, 1) - c0 - c1 - c2
            for x in np.roots([(c1 * np.conj(c12)).imag, (c0 * np.conj(c12) + c1 * np.conj(c2)).imag, (c0 * np.conj(c2)).imag]):
                if abs(x.imag) > 1e-9 * max(1, abs(x)): continue
                x = x.real; y = (-(c0 + c1 * x) / (c2 + c12 * x)).real; k = kp0.copy(); k[i] += x; k[j] += y
                if np.all((k >= KPMIN) & (k <= KPMAX)) and (best is None or np.hypot(x, y) < best[0]): best = (np.hypot(x, y), i, j, k)
        if best:
            row['remote_pair'] = f'{30+best[1]},{30+best[2]}'; row['remote_norm'] = best[0] / kp0[0]
            row['remote_two_Kp_beyond_guard'] = count(Model(np.r_[P0[:10], best[3], P0[20:30]], tau), -0.05)[0] - 1
        else: row['remote_pair'] = 'no admissible real solution'
    rows.append(row); print(row, flush=True); pd.DataFrame(rows).to_csv('T12b_single_site_delay.csv', index=False)
