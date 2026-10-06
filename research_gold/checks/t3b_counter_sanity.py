import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
import numpy as np
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
def phase(m, s): return np.linalg.slogdet(m.delta(s))[0]
def count(m, re0, re1, im0, im1, n=600):
    c = [complex(re0, im0), complex(re1, im0), complex(re1, im1), complex(re0, im1), complex(re0, im0)]
    tot = 0.0
    for a, b in zip(c[:-1], c[1:]):
        ts = list(np.linspace(0, 1, n)); v = [phase(m, a + (b - a) * t) for t in ts]; k = 0
        while k < len(ts) - 1:
            d = np.angle(v[k + 1] / v[k])
            if abs(d) > 0.3 and ts[k + 1] - ts[k] > 1e-10:
                tm = .5 * (ts[k] + ts[k + 1]); ts.insert(k + 1, tm); v.insert(k + 1, phase(m, a + (b - a) * tm))
            else: tot += d; k += 1
    return tot / 2 / np.pi
m0 = Model(P0, np.zeros(10)); ev = la.eigvals(m0.A)
print('tau=0: eig Re>0.01:', int(np.sum(ev.real > 0.01)), ' max Re:', np.sort(ev.real)[-4:], '| counter:', round(count(m0, 0.01, 80, -600, 600), 3))
print('tau=0: eig in Re>-0.05:', int(np.sum(ev.real > -0.05)), '| counter:', round(count(m0, -0.05, 80, -600, 600), 3))
m = Model(P0, TAU)
for lo, hi in [(0.5, 10), (10, 60), (60, 200), (200, 600)]:
    print('40 ms, Re in (0.01,80), Im in', (lo, hi), '->', round(count(m, 0.01, 80, lo, hi), 3))
print('40 ms near real axis Im in (-0.5,0.5):', round(count(m, 0.01, 80, -0.5, 0.5), 3))
