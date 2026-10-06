"""T12: cooperative delay compensation with ALL K_I fixed (proposal reviewed 5 Oct): keep a chosen pole at uniform delay 44 ms by changing
only proportional gains.  (a) two sites: exact finite solution from det = c0 + c1 x + c2 y + c12 xy  -> real quadratic in x;
(b) all ten K_p: minimum-norm solution of the two real pole conditions (Newton on the exact multi-affine determinant).
IEEE-39 design A, exported 204-state DDE.  Spectral screen; events are run separately."""
import os, sys, pathlib, itertools
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
import model as MM
from model import Model, P0, TAU
from t1_spillover_ieee39 import roots
from t3_full_spectrum_count import count
H = 0.004; TAUN = TAU + H; KPMIN, KPMAX = 7.853981633974483, 125.66370614359172; GUARD = -0.05
base = Model(P0, TAU)
def detK(lam, kp):
    rb = la.solve(lam * base.eye - base.A0, MM.M['Bp'] * kp + MM.M['Bi'] * P0[20:30], check_finite=False)
    return np.linalg.det(np.eye(10) - (base.C @ rb) * np.exp(-lam * TAUN)[None, :])
def tracked_worst(p):
    m = Model(p, TAUN); w = -9
    for z in roots:
        try: w = max(w, m.refine(z)[0].real)
        except Exception: pass
    return w
fast = max([z for z in roots if z.imag > 25], key=lambda z: z.real); mid = min(roots, key=lambda z: abs(z.imag - 6.544))
kp0 = P0[10:20].copy(); cands = []
for tname, lam in [('1.04 Hz', mid), ('4.91 Hz', fast)]:
    scale = abs(detK(lam, kp0))
    for i, j in itertools.combinations(range(10), 2):
        def f(x, y):
            k = kp0.copy(); k[i] += x; k[j] += y; return detK(lam, k) / scale
        c0 = f(0, 0); c1 = f(1, 0) - c0; c2 = f(0, 1) - c0; c12 = f(1, 1) - c0 - c1 - c2
        A = (c1 * np.conj(c12)).imag; B = (c0 * np.conj(c12) + c1 * np.conj(c2)).imag; C = (c0 * np.conj(c2)).imag
        for x in np.roots([A, B, C]):
            if abs(x.imag) > 1e-9 * max(1, abs(x)): continue
            x = x.real; den = c2 + c12 * x
            if abs(den) < 1e-14: continue
            y = (-(c0 + c1 * x) / den).real; k = kp0.copy(); k[i] += x; k[j] += y
            res = abs(detK(lam, k)) / scale
            cands.append(dict(target=tname, kind='two Kp', sites=f'{30+i},{30+j}', dkp_i=x, dkp_j=y, in_bounds=bool(np.all((k >= KPMIN) & (k <= KPMAX))),
                              residual=res, norm_change=float(np.hypot(x, y) / kp0[0]), kp=k))
    # (b) all ten Kp, minimum norm
    x = np.zeros(10)
    for it in range(40):
        f0 = detK(lam, kp0 + x) / scale
        J = np.array([detK(lam, kp0 + x + np.eye(10)[q]) / scale - f0 for q in range(10)]); Jr = np.vstack([J.real, J.imag])
        step = Jr.T @ la.solve(Jr @ Jr.T, np.array([f0.real, f0.imag])); x = x - step
        if la.norm(step) < 1e-12: break
    k = kp0 + x
    cands.append(dict(target=tname, kind='ten Kp min-norm', sites='all', dkp_i=float(x.min()), dkp_j=float(x.max()), in_bounds=bool(np.all((k >= KPMIN) & (k <= KPMAX))),
                      residual=abs(detK(lam, k)) / scale, norm_change=float(la.norm(x) / kp0[0]), kp=k))
df = pd.DataFrame(cands); print('candidates', len(df), '| in bounds', int(df.in_bounds.sum()), '| by target', df[df.in_bounds].groupby(['target', 'kind']).size().to_dict(), flush=True)
ok = df[df.in_bounds & (df.residual < 1e-8)].copy()
ok['tracked_worst'] = [tracked_worst(np.r_[P0[:10], k, P0[20:30]]) for k in ok.kp]
ok = ok.sort_values('tracked_worst'); print(ok.drop(columns='kp').round(4).head(14).to_string(index=False), flush=True)
rows = []
sel = pd.concat([ok[ok.kind == 'ten Kp min-norm'], ok[(ok.kind == 'two Kp') & (ok.tracked_worst < GUARD)].sort_values('norm_change').head(5)])
for _, r in sel.iterrows():
    p = np.r_[P0[:10], r.kp, P0[20:30]]; m = Model(p, TAUN); nu = count(m, 0.01)[0]; ng = count(m, GUARD)[0] - 1
    rows.append(dict(target=r.target, kind=r.kind, sites=r.sites, unstable=nu, beyond_guard=ng, norm_change=r.norm_change, tracked_worst=r.tracked_worst,
                     **{f'kp{30+q}': r.kp[q] for q in range(10)}))
    print({k: (round(float(v), 4) if not isinstance(v, str) else v) for k, v in rows[-1].items() if not k.startswith('kp3')}, flush=True)
    pd.DataFrame(rows).to_csv('T12_kp_only_cooperative.csv', index=False)
ok.drop(columns='kp').to_csv('T12_all_inbound_candidates.csv', index=False)
