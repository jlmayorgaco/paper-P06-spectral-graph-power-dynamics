"""T6: the interpolation node is a free design choice. Protect the heavily damped 1.04 Hz mode instead of the 4.91 Hz one:
same full-spectrum stability? how much less K_I is spent (K_I sets the low-frequency moments)?"""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
from t1_spillover_ieee39 import transport, roots
from t3_full_spectrum_count import count
mid = min(roots, key=lambda z: abs(z.imag - 6.544)); print('node', mid)
rows = []
for hms in [2, 4, 6, 8]:
    h = hms * 1e-3; p = P0.copy()
    for i in range(10): p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], mid, h)
    m = Model(p, TAU + h); nu, _ = count(m, 0.01); ng, _ = count(m, -0.05)
    rows.append(dict(node_Hz=mid.imag / 2 / np.pi, tau_ms=40 + hms, unstable_roots=nu, right_of_guard_excl_gauge=ng - 1, Kp=p[10], Ki=p[20]))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv('T6_node_choice.csv', index=False)
a, w = mid.real, mid.imag
print('budget h* [ms] for this node:', 1e3 * np.arctan2(w * P0[20], abs(mid) ** 2 * P0[10] + a * P0[20]) / w)
