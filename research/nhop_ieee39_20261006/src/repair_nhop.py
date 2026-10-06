"""E3/E4 engine: generalisation of research/feedback_cycle_fullmodel_20261005/src/repair.py (engine v2) to matrix gains.
Same logic: LP (min worst-root t) -> SLSQP min-norm step with target, trust region (cap 0.1, halve on rejection), exact exponential-delay
catalogue recomputed after every trial step, acceptance iff worst catalogued root improves by > 1e-6, best iterate kept, <=60 iterations,
stop when worst Re <= -SIGMA-1e-4 or radius < 1e-3. Success is decided by the full-spectrum counter (N_margin = 0, gauge excluded).
Variables: log Kp_ii, log KI_ii in [0.25,4]x nominal; u_ij, v_ij in [-1,1] (kp_ij = KP0*u_ij, kI_ij = KI0*v_ij), j in N_i^(n), j != i."""
import sys, pathlib, time, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from nhop import *
from scipy.optimize import linprog, minimize
STEP = 0.1; BUF = 0.003; MAXIT = 60
BP, BI = M['Bp'], M['Bi']


def active_roots(m, floor=-1.5):
    cat, _ = catalog(m, re_floor=-6, fmax=20)
    return [z for z in cat if z.imag > 0.3 and z.real > floor]


def build_gains(x, off):
    """x = [logKp_ii(10), logKI_ii(10), u(offs), v(offs)]"""
    Kp = np.diag(np.exp(x[:10])); KI = np.diag(np.exp(x[10:20])); no = len(off)
    for k, (i, j) in enumerate(off):
        Kp[i, j] = KP0 * x[20 + k]; KI[i, j] = KI0 * x[20 + no + k]
    return Kp, KI


def sens_all(m, lam, x, off):
    """exact simple-root sensitivities d lam / d(variable) for all variables (same formula as gain_authority.sens, matrix gains)."""
    r, l, _ = null_vectors(m, lam); e = np.exp(-lam * m.tau); Ce = m.C @ r
    Ds_r = r + m.B @ ((m.tau * e) * Ce)                      # Delta'(lam) r
    den = np.vdot(l, Ds_r); w = e * Ce / den
    ap = np.array([np.vdot(l, BP[:, i]) for i in range(10)]); ai = np.array([np.vdot(l, BI[:, i]) for i in range(10)])
    Sp = np.outer(ap, w); Si = np.outer(ai, w)               # d lam / d Kp_ij , d lam / d KI_ij
    kd = np.exp(x[:20])
    return np.r_[[Sp[i, i] * kd[i] for i in range(10)], [Si[i, i] * kd[10 + i] for i in range(10)],
                 [Sp[i, j] * KP0 for i, j in off], [Si[i, j] * KI0 for i, j in off]]


def fd_check(m0, lam, x, off, tau, rho, h=1e-5):
    g = sens_all(m0, lam, x, off); out = {}
    cases = {'diag_logKp0': 0, 'diag_logKI0': 10}
    if off: cases['offdiag_u0'] = 20; cases['offdiag_v0'] = 20 + len(off)
    for name, k in cases.items():
        zs = []
        for sgn in (1, -1):
            xx = x.copy(); xx[k] += sgn * h; Kp, KI = build_gains(xx, off)
            zs.append(DModel(Kp, KI, tau, rho).refine(lam)[0])
        fd = (zs[0] - zs[1]) / (2 * h); out[name] = float(abs(fd - g[k]) / abs(fd))
    return out


def baseline_x(off): return np.r_[np.log(KPB) * np.ones(10), np.log(KIB) * np.ones(10), np.zeros(2 * len(off))]


def repair_nhop(dist, n, tau=0.044, rho=0.875, maxit=MAXIT, label='', verbose=True):
    t0 = time.time(); tauv = np.full(10, tau)
    off = [(i, j) for i in range(10) for j in range(10) if j != i and dist[i, j] <= n]; no = len(off)
    x = baseline_x(off); nv = 20 + 2 * no
    lo = np.r_[np.log(.25 * KP0) * np.ones(10), np.log(.25 * KI0) * np.ones(10), -np.ones(2 * no)]
    hi = np.r_[np.log(4 * KP0) * np.ones(10), np.log(4 * KI0) * np.ones(10), np.ones(2 * no)]
    Kp, KI = build_gains(x, off); m = DModel(Kp, KI, tauv, rho); act = active_roots(m)
    worst = max(act, key=lambda z: z.real); cur_w = worst.real; radius = STEP; best_x, best_w, best_z = x.copy(), np.inf, None; hist = []; fd = None
    fd = fd_check(m, m.refine(worst)[0], x, off, tauv, rho)
    for it in range(maxit + 1):
        nb = int(sum(z.real > -SIGMA for z in act))
        hist.append(dict(it=it, worst_re=cur_w, worst_f_hz=worst.imag / 2 / np.pi, n_active_beyond=nb, radius=radius, elapsed=time.time() - t0))
        if verbose: print(f'[{label}] it {it}: worst {cur_w:+.4f} @ {worst.imag/2/np.pi:.3f} Hz, beyond {nb}, radius {radius:.4f}', flush=True)
        if cur_w < best_w: best_w, best_x, best_z = cur_w, x.copy(), worst
        if cur_w <= -SIGMA - 1e-4 or it == maxit or radius < 1e-3: break
        G = np.array([sens_all(m, z, x, off).real for z in act]); R0 = np.array([z.real for z in act])
        lb = np.maximum(-radius, lo - x); ub = np.minimum(radius, hi - x)
        res = linprog(np.r_[np.zeros(nv), 1.], A_ub=np.c_[G, -np.ones(len(R0))], b_ub=-R0, bounds=list(zip(lb, ub)) + [(None, None)], method='highs')
        if res.status != 0:
            hist[-1]['note'] = 'LP infeasible'; break
        target = max(-SIGMA - BUF, res.x[-1] + 1e-6)
        cons = {'type': 'ineq', 'fun': lambda d: target - (R0 + G @ d), 'jac': lambda d: -G}
        q = minimize(lambda d: .5 * d @ d, res.x[:nv], jac=lambda d: d, bounds=list(zip(lb, ub)), constraints=[cons], method='SLSQP', options=dict(maxiter=200, ftol=1e-12))
        d = q.x if q.success or np.all(target - (R0 + G @ q.x) > -1e-6) else res.x[:nv]
        trial = np.clip(x + d, lo, hi); Kt, It = build_gains(trial, off); mt = DModel(Kt, It, tauv, rho); at = active_roots(mt); wt = max(at, key=lambda z: z.real)
        if wt.real < cur_w - 1e-6:
            x, m, act, worst, cur_w = trial, mt, at, wt, wt.real; hist[-1]['accepted'] = True
        else:
            radius /= 2; hist[-1]['accepted'] = False
    Kp, KI = build_gains(best_x, off); mb = DModel(Kp, KI, tauv, rho); nu, nm, wind = spectrum_counts(mb)
    Kp0m = np.diag(KPB * np.ones(10)); KI0m = np.diag(KIB * np.ones(10))
    effort = float(la.norm(Kp - Kp0m) / KP0 + la.norm(KI - KI0m) / KI0)
    final = dict(n=n, tau_ms=tau * 1000, rho=rho, n_offdiag_vars=2 * no, N_unstable=nu, N_margin=nm, success=bool(nm == 0), iterations=len(hist) - 1,
                 worst_catalog_re=float(best_w), worst_catalog_f_hz=float(best_z.imag / 2 / np.pi), effort=effort, max_noninteger_winding=float(wind),
                 fd_rel_err=fd, runtime_s=time.time() - t0, stop_note=hist[-1].get('note', ''))
    return Kp, KI, hist, final


def write_toml(path, rho, tau, Kp, KI):
    def mat(M_): return '[\n' + ''.join('  [' + ', '.join(repr(float(v)) for v in row) + '],\n' for row in M_) + ']'
    pathlib.Path(path).write_text(f'rho = [{", ".join([repr(float(rho))] * 10)}]\ntau = {tau:.3f}\nKp = {mat(Kp)}\nKI = {mat(KI)}\n')
