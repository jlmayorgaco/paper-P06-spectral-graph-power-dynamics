"""Sparse PLL retuning engine (Phases I/J). Trust-region sensitivity step + EXACT full-model recomputation after every finite step.
Variables: log Kp_i and/or log Ki_i on a given support. Frozen: tau, sigma_req=0.05, gain bounds 0.25-4x nominal, per-iteration step cap 0.1 (log), <= 60 iterations.
min 1/2 |delta|^2  s.t.  Re(lam_r) + g_r . delta <= -(sigma_req + buffer) for every catalogued root with Re > -1.5,  bounds, step cap."""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
from gain_authority import sens
from scipy.optimize import linprog, minimize
from t3_full_spectrum_count import count
STEP = 0.1; BUF = 0.003; MAXIT = 60
def active_roots(m, floor=-1.5):
    cat, _ = catalog(m, re_floor=-6, fmax=20)
    return [z for z in cat if z.imag > 0.3 and z.real > floor]
def repair(sites, kind, tau_ms=44., p0=None, maxit=MAXIT, verbose=True, label='', tied=False):
    """sites: iterable of bus numbers (30..39). kind in {'ki','kp','joint'}."""
    p = (P0 if p0 is None else p0).copy(); tau = np.full(10, tau_ms / 1000); idx = [s - 30 for s in sites]
    var = []   # (tuple of param indices moved together, lo, hi)
    if tied:
        if kind in ('kp', 'joint'): var.append((tuple(10 + i for i in idx), KPLO, KPHI))
        if kind in ('ki', 'joint'): var.append((tuple(20 + i for i in idx), KILO, KIHI))
    else:
        for i in idx:
            if kind in ('kp', 'joint'): var.append(((10 + i,), KPLO, KPHI))
            if kind in ('ki', 'joint'): var.append(((20 + i,), KILO, KIHI))
    hist = []; radius = STEP; halv = 0; t0 = time.time(); best_p = None; best_w = np.inf
    m = Model(p, tau); act = active_roots(m); worst = max(act, key=lambda z: z.real); cur_w = worst.real
    for it in range(maxit + 1):
        hist.append(dict(label=label, kind=kind, sites=' '.join(map(str, sites)), iteration=it, worst_re=cur_w, worst_f_hz=worst.imag / 2 / np.pi, n_active_gt_margin=int(sum(z.real > -SIGMA for z in act)), radius=radius,
                         Kp=' '.join(f'{x:.3f}' for x in p[10:20]), Ki=' '.join(f'{x:.3f}' for x in p[20:30]), elapsed_s=time.time() - t0))
        if verbose: print(f'[{label}] it {it}: worst Re={cur_w:+.4f} @ {worst.imag / 2 / np.pi:.3f} Hz, beyond margin={hist[-1]["n_active_gt_margin"]}, radius={radius:.4f}', flush=True)
        if cur_w < best_w: best_w, best_p = cur_w, p.copy()
        if cur_w <= -SIGMA - 1e-4 or it == maxit or radius < 1e-3: break
        G = []; R0 = []
        for z in act:
            dKp, dKi = sens(m, z); G.append([sum((dKp[pi - 10] if pi < 20 else dKi[pi - 20]).real for pi in pis) for pis, lo, hi in var]); R0.append(z.real)
        G = np.array(G); R0 = np.array(R0); nv = len(var)
        lb = np.array([max([-radius] + [np.log(lo / p[pi]) for pi in pis]) for pis, lo, hi in var]); ub = np.array([min([radius] + [np.log(hi / p[pi]) for pi in pis]) for pis, lo, hi in var])
        res = linprog(np.r_[np.zeros(nv), 1.], A_ub=np.c_[G, -np.ones(len(R0))], b_ub=-R0, bounds=list(zip(lb, ub)) + [(None, None)], method='highs')
        if res.status != 0: hist[-1]['note'] = 'LP infeasible'; break
        tstar = res.x[-1]; target = max(-SIGMA - BUF, tstar + 1e-6)
        cons = {'type': 'ineq', 'fun': lambda d: target - (R0 + G @ d), 'jac': lambda d: -G}
        q = minimize(lambda d: .5 * d @ d, res.x[:nv], jac=lambda d: d, bounds=list(zip(lb, ub)), constraints=[cons], method='SLSQP', options=dict(maxiter=200, ftol=1e-12))
        d = q.x if q.success or np.all(target - (R0 + G @ q.x) > -1e-6) else res.x[:nv]
        hist[-1]['predicted_worst'] = float(np.max(R0 + G @ d)); hist[-1]['step_norm'] = float(np.linalg.norm(d))
        trial = p.copy()
        for (pis, lo, hi), di in zip(var, d):
            for pi in pis: trial[pi] = float(np.clip(trial[pi] * np.exp(di), lo, hi))
        mt = Model(trial, tau); at = active_roots(mt); wt = max(at, key=lambda z: z.real)
        if wt.real < cur_w - 1e-6: p, m, act, worst, cur_w = trial, mt, at, wt, wt.real; hist[-1]['accepted'] = True
        else: radius /= 2; hist[-1]['accepted'] = False; hist[-1]['trial_worst'] = wt.real
    p = best_p
    m = Model(p, tau); nu, e1 = count(m, .01); ng, e2 = count(m, -SIGMA)
    final = dict(label=label, kind=kind, sites=' '.join(map(str, sites)), N_unstable=nu, N_margin_excl_gauge=ng - 1, max_noninteger_winding=max(e1, e2), iterations=len(hist) - 1,
                 success=bool(ng - 1 == 0), Kp=' '.join(f'{x:.4f}' for x in p[10:20]), Ki=' '.join(f'{x:.4f}' for x in p[20:30]), worst_re_catalog=best_w)
    if verbose: print('FINAL', final, flush=True)
    return p, hist, final
