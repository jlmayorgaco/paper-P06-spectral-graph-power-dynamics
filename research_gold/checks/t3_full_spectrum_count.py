"""T3: banded argument-principle root count on det of the full 204x204 characteristic matrix (exact exponentials).
Region: Re in (re0, 80), |Im| < 600 rad/s (95 Hz), split in frequency bands (a single large rectangle aliases).
Floating point with adaptive bisection; validated against eig() at tau=0. NOT an interval certificate."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
from t1_spillover_ieee39 import transport, roots
BANDS = [(0.5, 10), (10, 60), (60, 200), (200, 600)]
def phase(m, s): return np.linalg.slogdet(m.delta(s))[0]
def box(m, re0, re1, im0, im1, n=600):
    c = [complex(re0, im0), complex(re1, im0), complex(re1, im1), complex(re0, im1), complex(re0, im0)]; tot = 0.0
    for a, b in zip(c[:-1], c[1:]):
        ts = list(np.linspace(0, 1, n)); v = [phase(m, a + (b - a) * t) for t in ts]; k = 0
        while k < len(ts) - 1:
            d = np.angle(v[k + 1] / v[k])
            if abs(d) > 0.3 and ts[k + 1] - ts[k] > 1e-10:
                tm = .5 * (ts[k] + ts[k + 1]); ts.insert(k + 1, tm); v.insert(k + 1, phase(m, a + (b - a) * tm))
            else: tot += d; k += 1
    return tot / 2 / np.pi
def count(m, re0):
    w = [box(m, re0, 80, -0.5, 0.5)] + [box(m, re0, 80, lo, hi) for lo, hi in BANDS]
    return int(round(w[0])) + 2 * sum(int(round(x)) for x in w[1:]), max(abs(x - round(x)) for x in w)
if __name__ == '__main__':
    m0 = Model(P0, np.zeros(10)); ev = la.eigvals(m0.A)
    print('VALIDATION tau=0  eig(Re>0.01)=', int(np.sum(ev.real > 0.01)), 'counter=', count(m0, 0.01),
          '| eig(Re>-0.05)=', int(np.sum(ev.real > -0.05)), 'counter=', count(m0, -0.05), flush=True)
    slow = [z for z in roots if z.imag < 1][0]; fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real)
    rows = []
    for name, lam, hs in [('none', None, [0, 1, 2, 4]), ('protect_slow', slow, [1, 2, 4]), ('protect_fast', fast, [2, 4, 6, 8, 8.9])]:
        for hms in hs:
            h = hms * 1e-3; p = P0.copy()
            if lam is not None:
                for i in range(10): p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], lam, h)
            m = Model(p, TAU + h); nu, e1 = count(m, 0.01); ng, e2 = count(m, -0.05)
            rows.append(dict(strategy=name, tau_ms=40 + hms, unstable_roots=nu, roots_right_of_guard_excl_gauge=ng - 1,
                             max_noninteger_winding=round(max(e1, e2), 6), Kp=p[10], Ki=p[20]))
            print({k: (round(float(v), 4) if not isinstance(v, str) else v) for k, v in rows[-1].items()}, flush=True)
    pd.DataFrame(rows).to_csv('T3_full_spectrum_count.csv', index=False)
