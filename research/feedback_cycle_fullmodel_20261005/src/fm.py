"""Campaign library: thin read-only wrapper on the frozen exact-delay model (experiments/interaction_decision_20261004/model.py)."""
import os, sys, pathlib
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd
from scipy import linalg as la
CAMP = pathlib.Path(__file__).resolve().parents[1]
ROOT = CAMP.parents[1]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
sys.path.insert(0, str(ROOT / 'research_gold' / 'checks'))
import model as MM
from model import Model, P0, TAU, M
KP0, KI0 = 2 * np.pi * 5, (2 * np.pi * 5) ** 2 / 4
KPB, KIB = P0[10], P0[20]
BUSES = np.arange(30, 40)
KPLO, KPHI, KILO, KIHI = .25 * KP0, 4 * KP0, .25 * KI0, 4 * KI0
SIGMA = 0.05
def design(rho=0.875, kp=None, ki=None):
    p = P0.copy(); p[:10] = rho
    if kp is not None: p[10:20] = kp
    if ki is not None: p[20:30] = ki
    return p
def pll_index():
    """indices (theta, omega, xi) of the ten PLLs in the 204-state vector, from the exported port table."""
    ports = MM.PORTS
    return [(int(r.pll_angle_index), int(r.pll_frequency_index), int(r.pll_frequency_index) + 1) for r in ports.itertuples()]
def pad_seeds(m, order=6):
    """Pade(order) augmented eigenvalues: SEEDS ONLY (never reported as roots)."""
    from scipy.signal import tf2ss
    from math import factorial
    c = np.array([factorial(2 * order - k) * factorial(order) / (factorial(2 * order) * factorial(k) * factorial(order - k)) for k in range(order + 1)])
    A_, B_, C_, D_ = tf2ss((c * (-1.) ** np.arange(order + 1))[::-1], c[::-1])
    n = 204; N = n + order * 10; A = np.zeros((N, N)); A[:n, :n] = m.A0 + m.B[:, :0].sum() * 0
    A[:n, :n] = m.A0 + D_[0, 0] * sum(np.outer(m.B[:, i], m.C[i]) for i in range(10))
    for i in range(10):
        if m.tau[i] < 1e-12: raise ValueError('positive delay needed')
        ss = slice(n + i * order, n + (i + 1) * order)
        A[:n, ss] = np.outer(m.B[:, i], C_.ravel()); A[ss, :n] = np.outer(B_.ravel(), m.C[i]) / m.tau[i]; A[ss, ss] = A_ / m.tau[i]
    return la.eigvals(A)
def catalog(m, re_floor=-60, fmax=60., order=6):
    """All exact-delay roots with Im>=0 found by refining Pade seeds with the exact characteristic. Not a completeness proof."""
    seeds = pad_seeds(m, order)
    seeds = seeds[(seeds.real > re_floor) & (seeds.imag >= -1e-6) & (seeds.imag < 2 * np.pi * fmax)]
    out = []; res = []
    for s in sorted(seeds, key=lambda z: -z.real):
        try: z, e = m.refine(s)
        except Exception: continue
        if abs(z) < 1e-5: continue
        if z.imag < -1e-6: continue
        if all(abs(z - q) > 1e-6 for q in out): out.append(z); res.append(e)
    idx = np.argsort([-q.real for q in out]); return np.array(out)[idx], np.array(res)[idx]
def null_vectors(m, lam):
    """right/left null vectors of the full 204x204 characteristic Delta(lam) (exact exponentials)."""
    D = m.delta(lam); U, S, Vh = la.svd(D); return Vh[-1].conj(), U[:, -1], S[-1]
def participation(m, lam, groups):
    r, l, _ = null_vectors(m, lam); p = np.abs(l.conj() * r); p /= p.sum()
    return {k: float(p[ix].sum()) for k, ix in groups.items()}
def groups_for(tag='0p875'):
    """0-based state groups from the Julia state map: SG, GFL, PLL (theta,omega,xi) and per-site sets."""
    sm = pd.read_csv(CAMP / 'raw' / 'jac' / f'statemap_rho{tag}.csv'); sm['i'] = sm['index'] - 1
    g = {'SG': sm[sm.kind == 'SG'].i.to_numpy(), 'GFL': sm[sm.kind == 'GFL'].i.to_numpy()}
    pll = [x for t in pll_index() for x in t]; g['PLL'] = np.array(pll)
    for k, b in enumerate(BUSES):
        g[f'PLL{b}'] = np.array(pll_index()[k]); g[f'site{b}'] = sm[sm.site == b].i.to_numpy()
    return g
