"""T1 infrastructure: distributed-gain exact-delay IEEE-39 characteristic, root refinement, full-spectrum counter,
communication graph G_c, open-PLL operator G(s).  Reuses (read-only) research/feedback_cycle_fullmodel_20261005/src/fm.py."""
import os, sys, pathlib
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = pathlib.Path(__file__).resolve().parent
CAMP = HERE.parent
ROOT = CAMP.parents[1]
sys.path.insert(0, str(ROOT / 'research' / 'feedback_cycle_fullmodel_20261005' / 'src'))
from fm import *          # np, pd, la, Model, P0, M, KP0, KI0, KPB, KIB, SIGMA, pll_index, catalog, ...
from scipy.sparse.csgraph import shortest_path, connected_components
CAMP = HERE.parent          # re-bind AFTER the fm import (fm defines its own CAMP)
TF = 1 / (600 * np.pi)
EYE = np.eye(204)


class DModel:
    """Distributed PLL gains. Kp, KI : 10x10 (row = receiving PLL, col = measured detector)."""
    def __init__(self, Kp, KI, tau=0.044, rho=None):
        p = P0.copy()
        if rho is not None: p[:10] = rho
        self.base = Model(p, np.full(10, tau) if np.isscalar(tau) else tau)
        self.A0, self.C, self.tau = self.base.A0, self.base.C, self.base.tau
        self.Kp, self.KI = np.array(Kp, float), np.array(KI, float)
        self.B = M['Bp'] @ self.Kp + M['Bi'] @ self.KI           # 204 x 10

    def delta(self, s):
        return s * EYE - self.A0 - self.B @ (np.exp(-s * self.tau)[:, None] * self.C)

    def reduced(self, s):
        """K(s) = I - E C (sI-A0)^-1 B  (10x10); det Delta = det(sI-A0) det K."""
        rb = la.solve(s * EYE - self.A0, self.B, check_finite=False)
        return np.eye(10) - np.exp(-s * self.tau)[:, None] * (self.C @ rb), rb

    def refine(self, s, tol=1e-10, maxit=40):
        """Newton on the smallest singular value (u^H K v) of the reduced characteristic; exact exponentials."""
        s = complex(s)
        for _ in range(maxit):
            lu = la.lu_factor(s * EYE - self.A0, check_finite=False)
            rb = la.lu_solve(lu, self.B, check_finite=False); e = np.exp(-s * self.tau)
            K = np.eye(10) - e[:, None] * (self.C @ rb)
            Ul, sv, vh = la.svd(K, check_finite=False); v = vh[-1].conj(); u = Ul[:, -1]
            Ks = e[:, None] * (self.C @ la.lu_solve(lu, rb, check_finite=False)) + (self.tau * e)[:, None] * (self.C @ rb)
            step = -np.vdot(u, K @ v) / np.vdot(u, Ks @ v)
            if abs(step) < 1e-10 * max(1, abs(s)):
                full = rb @ v; Dl = self.delta(s); res = la.norm(Dl @ full) / (la.norm(Dl) * la.norm(full))
                if res < tol: return s, float(res)
            if abs(step) > 1: step /= abs(step)
            s += step
        raise RuntimeError(f'refine failed s={s}')

    def null_right(self, s):
        """right null vector of Delta(s) (204, unit 2-norm): v = (sI-A0)^-1 B w, K w = 0."""
        K, rb = self.reduced(s); w = la.svd(K)[2][-1].conj(); v = rb @ w; return v / la.norm(v)


def diag_model(tau, rho=None, kp=None, ki=None):
    kp = KPB * np.ones(10) if kp is None else kp; ki = KIB * np.ones(10) if ki is None else ki
    return DModel(np.diag(kp), np.diag(ki), tau, rho)


# ---------------- full-spectrum counter (generalised from research_gold/checks/t3_full_spectrum_count.py; same bands)
BANDS = [(0.5, 10), (10, 60), (60, 200), (200, 600)]
def _phase(m, s): return np.linalg.slogdet(m.delta(s))[0]
def _box(m, re0, re1, im0, im1, n=600):
    c = [complex(re0, im0), complex(re1, im0), complex(re1, im1), complex(re0, im1), complex(re0, im0)]; tot = 0.0
    for a, b in zip(c[:-1], c[1:]):
        ts = list(np.linspace(0, 1, n)); v = [_phase(m, a + (b - a) * t) for t in ts]; k = 0
        while k < len(ts) - 1:
            d = np.angle(v[k + 1] / v[k])
            if abs(d) > 0.3 and ts[k + 1] - ts[k] > 1e-10:
                tm = .5 * (ts[k] + ts[k + 1]); ts.insert(k + 1, tm); v.insert(k + 1, _phase(m, a + (b - a) * tm))
            else: tot += d; k += 1
    return tot / 2 / np.pi
def count(m, re0):
    """roots with Re > re0 (Re<80, |Im|<600): strip |Im|<0.5 once + 4 bands counted twice (conjugates)."""
    w = [_box(m, re0, 80, -0.5, 0.5)] + [_box(m, re0, 80, lo, hi) for lo, hi in BANDS]
    return int(round(w[0])) + 2 * sum(int(round(x)) for x in w[1:]), max(abs(x - round(x)) for x in w)
def spectrum_counts(m):
    """(N_unstable Re>0.01, N_margin Re>-0.05 excluding the gauge root at 0, max noninteger winding)"""
    nu, e1 = count(m, 0.01); ng, e2 = count(m, -SIGMA)
    return nu, ng - 1, max(e1, e2)


# ---------------- open-PLL operator
PLL = pll_index(); TH = np.array([t[0] for t in PLL]); OM = np.array([t[1] for t in PLL]); XI = np.array([t[2] for t in PLL])
PLLALL = np.array([x for t in PLL for x in t]); REST = np.setdiff1d(np.arange(204), PLLALL)
class OpenPLL:
    def __init__(self, rho=None):
        mm = Model(P0 if rho is None else np.r_[np.full(10, rho), P0[10:]], np.full(10, .044)); self.A0, self.C = mm.A0, mm.C
        self.Arr = self.A0[np.ix_(REST, REST)]; self.Art = self.A0[np.ix_(REST, TH)]
        self.Ct = self.C[:, TH]; self.Cr = self.C[:, REST]
        self.chk = dict(A_rest_rows_from_omega_xi=float(abs(self.A0[np.ix_(REST, np.r_[OM, XI])]).max()),
                        C_cols_omega_xi=float(abs(self.C[:, np.r_[OM, XI]]).max()),
                        A_PLLrows_from_rest=float(abs(self.A0[np.ix_(PLLALL, REST)]).max()))
        self.nr = len(REST)
    def G(self, s): return self.Ct + self.Cr @ la.solve(s * np.eye(self.nr) - self.Arr, self.Art, check_finite=False)


# ---------------- communication graph
def port_admittance():
    Y = pd.read_csv(ROOT / 'experiments/graph_gsp_codesign_20261003/model/Y.csv').to_numpy()
    Yc = np.zeros((10, 10), complex); viol = 0.0
    for i in range(10):
        for j in range(10):
            B = Y[2*i:2*i+2, 2*j:2*j+2]
            viol = max(viol, abs(B[0, 0] - B[1, 1]), abs(B[0, 1] + B[1, 0]))   # form [[a,-b],[b,a]]
            Yc[i, j] = B[0, 0] + 1j * B[1, 0]
    info = dict(block_form_violation_max=float(viol), convention='2x2 block [[a,-b],[b,a]] -> a+jb (b = lower-left entry)',
                Yc_symmetry_err=float(abs(Yc - Yc.T).max()), rcond=float(1 / np.linalg.cond(Yc)))
    return Yc, info

def comm_graph(bridge=False):
    Yc, info = port_admittance()
    if info['rcond'] > 1e-12: Z = np.linalg.inv(Yc); info['Z'] = 'inverse'
    else: Z = np.linalg.pinv(Yc); info['Z'] = 'pinv'
    R = np.abs(np.diag(Z)[:, None] + np.diag(Z)[None, :] - 2 * Z)
    A = np.zeros((10, 10), int)
    for i in range(10):
        for j in [j for j in np.argsort(R[i], kind='stable') if j != i][:2]: A[i, j] = A[j, i] = 1
    info['bridge_added'] = None
    if bridge:   # AMENDMENT_01 sensitivity variant only: join components by the single minimum-r cross pair
        ncomp, lab = connected_components(A)
        while ncomp > 1:
            cand = [(R[i, j], i, j) for i in range(10) for j in range(10) if lab[i] != lab[j]]; _, i, j = min(cand)
            A[i, j] = A[j, i] = 1; info['bridge_added'] = (int(i), int(j)); ncomp, lab = connected_components(A)
    ncomp, lab = connected_components(A); info['components'] = int(ncomp); info['labels'] = lab.tolist()
    D = shortest_path(A, unweighted=True)
    return dict(A=A, L=np.diag(A.sum(1)) - A, R=R, Z=Z, dist=D, diam=int(D[np.isfinite(D)].max()), info=info)  # diam = largest FINITE hop distance

def neighbourhoods(dist, n): return [np.where(dist[i] <= n)[0] for i in range(10)]   # inf distance never <= n
