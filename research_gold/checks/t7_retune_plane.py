"""T7: complete parametrisation of local PI retunes under a delay increment.
Any differentiable real retune is  k_p' = a k_p + k_I,  k_I' = -b k_p  with (a,b) in R^2; the loop-factor change is
-(s^2 - a s + b) k_p e^{-s tau} dh, so each mode moves by  d mu/dh = -(mu^2 - a mu + b) r_mu,  r_mu = sum_i k_p,i d mu/d k_I,i  (affine in (a,b)).
=> first-order safe set = intersection of half-planes (one per mode) in the (a,b) plane.
Special points: (0,0) moment/first-order compensation; (-k_I/k_p, 0) no retune; (2 Re lam, |lam|^2) node at lam."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy import linalg as la
from scipy.optimize import linprog
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
import model as MM
from model import Model, P0, TAU
from t1_spillover_ieee39 import roots
from t3_full_spectrum_count import count
base = Model(P0, TAU); H = 0.004; GUARD = -0.05
def muKi(mu):
    D = base.delta(mu); U, S, Vh = la.svd(D); v = Vh[-1].conj(); l = U[:, -1]
    den = np.vdot(l, base.derivative(mu) @ v)
    return np.array([np.vdot(l, MM.M['Bi'][:, i]) * np.exp(-mu * TAU[i]) * (base.C[i] @ v) / den for i in range(10)])
r = {mu: np.sum(P0[10:20] * muKi(mu)) for mu in roots}
# Re mu(h) ~ Re mu0 + H * [ -Re(mu^2 r) + a Re(mu r) - b Re(r) ]  <= GUARD
A = np.array([[H * (mu * r[mu]).real, -H * r[mu].real] for mu in roots]); c0 = np.array([mu.real - H * (mu ** 2 * r[mu]).real for mu in roots])
def worst(a, b): return float(np.max(c0 + A @ np.array([a, b])))
kp, ki = P0[10], P0[20]; fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real); mid = min(roots, key=lambda z: abs(z.imag - 6.544))
pts = {'moment (0,0)': (0.0, 0.0), 'no retune': (-ki / kp, 0.0), 'node 4.91 Hz': (2 * fast.real, abs(fast) ** 2), 'node 1.04 Hz': (2 * mid.real, abs(mid) ** 2)}
for k, (a, b) in pts.items(): print(f'{k:16s} a={a:9.3f} b={b:9.2f}  predicted worst Re = {worst(a, b):+.4f}')
# LP 1: best first-order margin.  LP 2: smallest |b| (least K_I drift, i.e. least change of the low-frequency moments) that keeps the guard.
res1 = linprog([0, 0, 1], A_ub=np.c_[A, -np.ones(len(roots))], b_ub=-c0, bounds=[(-200, 200), (-2000, 2000), (None, None)])
print('LP1 best margin: a,b,t =', res1.x)
res2 = linprog([0, 0, 1], A_ub=np.r_[np.c_[A, np.zeros(len(roots))], [[0, 1, -1], [0, -1, -1]]], b_ub=np.r_[GUARD - c0, 0, 0], bounds=[(-200, 200), (-2000, 2000), (0, None)])
print('LP2 least |b| with guard: a,b,|b| =', res2.x, res2.status)
rows = []
for name, (a, b) in {**pts, 'LP2 least-KI-drift': tuple(res2.x[:2])}.items():
    E = la.expm(H * np.array([[a, 1.0], [-b, 0.0]])); p = P0.copy()
    for i in range(10): p[10 + i], p[20 + i] = E @ np.array([P0[10 + i], P0[20 + i]])
    m = Model(p, TAU + H); nu, _ = count(m, 0.01); ng, _ = count(m, GUARD)
    rows.append(dict(retune=name, a=a, b=b, predicted_worst=worst(a, b), unstable_roots_44ms=nu, right_of_guard_excl_gauge=ng - 1, Kp=p[10], Ki=p[20]))
    print({k: (round(float(v), 4) if not isinstance(v, str) else v) for k, v in rows[-1].items()}, flush=True)
pd.DataFrame(rows).to_csv('T7_retune_plane.csv', index=False)
np.savez('T7_halfplanes.npz', A=A, c0=c0, roots=np.array(roots), H=H)
