"""E5: apply exact E1 assignments (eps<=1e-8), no corrector, count roots beyond the margin with the full counter.
variants: 'absolute' = row gains k_i = C_i^+ b_i (literal min-norm solution of the row conditions, baseline diagonal NOT kept);
          'delta'    = baseline diagonal + minimum-norm CHANGE dk_i = C_i^+ (b_i - C_i k0_i) (same residual, smallest departure from baseline)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
os_env = __import__('os').environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from e1_e2 import *
from multiprocessing import Pool


def build(D, mm, dist, n, variant):
    nb = neighbourhoods(dist, n); Kp = np.zeros((10, 10)); KI = np.zeros((10, 10)); worst = 0
    for i in range(10):
        C, b = rowsys(D, mm, i, nb[i]); k = len(nb[i])
        k0 = np.zeros(2 * k)
        if variant == 'delta': k0[list(nb[i]).index(i)] = KPB; k0[k + list(nb[i]).index(i)] = KIB
        x, res, rk = lsq(C, b - C @ k0); x = x + k0; worst = max(worst, np.linalg.norm(b - C @ x) / np.linalg.norm(b))
        Kp[i, nb[i]] = x[:k]; KI[i, nb[i]] = x[k:]
    return Kp, KI, worst


def job(args):
    gname, mm, n, variant, Kp, KI, lam_list, eps = args
    dm = DModel(Kp, KI, TAU44); nu, nm, w = spectrum_counts(dm)
    # target check: smallest singular value of the reduced characteristic at each lam*, and refined nearby root
    sv = []; dev = []
    for lam in lam_list:
        K, _ = dm.reduced(lam); sv.append(float(la.svdvals(K)[-1]))
        try: z, _ = dm.refine(lam); dev.append(abs(z - lam))
        except Exception: dev.append(np.nan)
    return dict(graph=gname, m=mm, n=n, variant=variant, eps_E1=eps, N_unstable=nu, N_margin=nm, max_noninteger_winding=float(w), spill_beyond_margin=bool(nm > 0),
                target_sigma_min_max=max(sv), target_refine_dev_max=float(np.nanmax(dev)), kp_diag_mean=float(np.mean(np.diag(Kp))), kI_diag_mean=float(np.mean(np.diag(KI))),
                kp_offdiag_max=float(np.max(np.abs(Kp - np.diag(np.diag(Kp))))), kI_offdiag_max=float(np.max(np.abs(KI - np.diag(np.diag(KI))))),
                kp_max=float(np.abs(Kp).max()), kI_max=float(np.abs(KI).max()))


if __name__ == '__main__':
    m, tg = targets(5); D = modal_data(m, tg); jobs = []
    for gname, br in (('frozen', False), ('bridge', True)):
        g = comm_graph(br)
        for mm in range(1, 6):
            for n in range(0, g['diam'] + 1):
                eps = free_nhop(D, mm, g['dist'], n)[0].max()
                if eps > TOL: continue
                for var in ('absolute', 'delta'):
                    Kp, KI, wres = build(D, mm, g['dist'], n, var); assert wres < 1e-8, (gname, mm, n, var, wres)
                    jobs.append((gname, mm, n, var, Kp, KI, [d['lam'] for d in D[:mm]], eps))
    print(len(jobs), 'jobs', flush=True)
    with Pool(min(20, len(jobs))) as P: out = P.map(job, jobs, chunksize=1)
    df = pd.DataFrame(out); df.to_csv(CAMP / 'derived' / 'E5_spillover.csv', index=False); print(df.to_string())
