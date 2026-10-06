"""T8: (a) first-order Taylor of the modal transport as a comparator (k_p + h(2a k_p + k_I), k_I - h|lam|^2 k_p) vs the exact transport
and vs the low-frequency rule (k_p + h k_I);  (b) the node-selection rule fixed on design A, applied blind to a second base design B.
RULE (fixed before looking at design B): among catalogued oscillatory modes (Im > 0.2 rad/s) whose toll-predicted worst root after +4 ms
is below the -0.05 guard, protect the one that leaves the largest K_I."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
import model as MM
from model import Model, P0, TAU, ROOTS
from t1_spillover_ieee39 import transport
from t3_full_spectrum_count import count
GUARD = -0.05
def catalog(m):
    out = []
    for s in ROOTS:
        try:
            r, _ = m.refine(s)
            if r.imag > 1e-6 and all(abs(r - q) > 1e-6 for q in out): out.append(r)
        except Exception: pass
    return sorted(out, key=lambda z: z.imag)
def gains(p0, rule, lam, h):
    p = p0.copy()
    for i in range(10):
        kp, ki = p0[10 + i], p0[20 + i]
        if rule == 'exact': p[10 + i], p[20 + i] = transport(kp, ki, lam, h)
        elif rule == 'taylor': p[10 + i], p[20 + i] = kp + h * (2 * lam.real * kp + ki), ki - h * abs(lam) ** 2 * kp
        elif rule == 'lowfreq': p[10 + i] = kp + h * ki
    return p
def muKi(m, p0, mu):
    D = m.delta(mu); U, S, Vh = la.svd(D); v = Vh[-1].conj(); l = U[:, -1]; den = np.vdot(l, m.derivative(mu) @ v)
    return np.array([np.vdot(l, MM.M['Bi'][:, i]) * np.exp(-mu * TAU[i]) * (m.C[i] @ v) / den for i in range(10)])
def select_node(p0, H=0.004):
    m = Model(p0, TAU); rts = catalog(m); sens = {mu: np.sum(p0[10:20] * muKi(m, p0, mu)) for mu in rts}; cand = []
    for lam in rts:
        if lam.imag < 0.2: continue
        pred = max(mu.real + H * (-(mu - lam) * (mu - lam.conjugate()) * sens[mu]).real for mu in rts if abs(mu - lam) > 1e-9)
        cand.append((lam, pred, transport(p0[10], p0[20], lam, H)[1]))
    safe = [c for c in cand if c[1] <= GUARD]
    return (max(safe, key=lambda c: c[2]) if safe else None), cand, rts
rows = []
def run(design, p0, name, rule, lam, taus):
    for t in taus:
        h = (t - 40) * 1e-3; p = gains(p0, rule, lam, h); m = Model(p, TAU + h); nu, _ = count(m, 0.01); ng, _ = count(m, GUARD)
        rows.append(dict(design=design, strategy=name, tau_ms=t, unstable_roots=nu, beyond_guard_excl_gauge=ng - 1, Kp=p[10], Ki=p[20]))
        print(rows[-1], flush=True); pd.DataFrame(rows).to_csv('T8_taylor_and_second_design.csv', index=False)
if __name__ == '__main__':
    # ---- (a) design A comparators
    selA, candA, rtsA = select_node(P0); print('RULE on design A selects', selA, flush=True)
    fast = max([z for z in rtsA if z.imag > 25], key=lambda z: z.real); mid = min(rtsA, key=lambda z: abs(z.imag - 6.544))
    run('A', P0, 'taylor node 4.91 Hz', 'taylor', fast, [42, 44, 46, 48]); run('A', P0, 'taylor node 1.04 Hz', 'taylor', mid, [42, 44, 46, 48])
    run('A', P0, 'low-frequency rule', 'lowfreq', None, [46])
    # ---- (b) design B: rho = 0.85 at every site, baseline gains
    PB = P0.copy(); PB[:10] = 0.85
    print('design B base count (unstable, beyond guard):', count(Model(PB, TAU), 0.01)[0], count(Model(PB, TAU), GUARD)[0] - 1, flush=True)
    selB, candB, rtsB = select_node(PB); print('catalog B', [f'{z.real:+.3f}{z.imag:+.2f}j' for z in rtsB]); print('RULE on design B selects', selB, flush=True)
    pd.DataFrame([dict(node=str(c[0]), predicted_worst=c[1], Ki_after_4ms=c[2]) for c in candB]).to_csv('T8_designB_candidates.csv', index=False)
    run('B', PB, 'no retune', 'none', None, [42, 44]); run('B', PB, 'low-frequency rule', 'lowfreq', None, [42, 44])
    if selB: run('B', PB, f'exact, rule-selected node {selB[0].imag/2/np.pi:.2f} Hz', 'exact', selB[0], [42, 44, 46])
    fastB = max([z for z in rtsB if z.imag > 25], key=lambda z: z.real); run('B', PB, f'exact, least-damped cluster node {fastB.imag/2/np.pi:.2f} Hz', 'exact', fastB, [44])
