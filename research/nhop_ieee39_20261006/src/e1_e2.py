"""E1 (free n-hop modal authority) and E2 (architecture comparison). Graphs: frozen G_c (primary) and AMENDMENT_01 bridge (sensitivity)."""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from nhop import *
TAU44 = .044; LAMRE = -0.30; TOL = 1e-8


def targets(mmax=5):
    """m rightmost PLL-family roots (PLL participation >= 0.5 of the Delta null vector) at 44 ms, baseline gains."""
    m = diag_model(TAU44); ref = Model(P0, np.full(10, TAU44)); cat, _ = catalog(ref, re_floor=-3, fmax=20); g = groups_for(); out = []
    for z in cat:
        if z.imag < 0.3: continue
        z, _ = m.refine(z); part = participation(ref, z, {'PLL': g['PLL']})['PLL']
        if part >= 0.5: out.append((z, part))
        if len(out) == mmax: break
    return m, out


def modal_data(m, tgts):
    """per target: lam*, q (theta comps of right null vector, unit 2-norm), y = G(lam*) q, z_j = exp(-lam* tau_j) y_j, b_i."""
    op = OpenPLL(); D = []
    for z, part in tgts:
        v = m.null_right(z); q = v[TH]; q = q / la.norm(q)
        lam = LAMRE + 1j * z.imag; y = op.G(lam) @ q; zz = np.exp(-lam * m.tau) * y; b = lam**2 * (1 + TF * lam) * q
        D.append(dict(root=z, part=part, lam=lam, q=q, y=y, z=zz, b=b))
    return D


def realify(Cc, bc):
    return np.vstack([Cc.real, Cc.imag]), np.r_[bc.real, bc.imag]


def rowsys(D, mm, i, nbr):
    """real system C_i k_i = b_i for site i, targets 0..mm-1, unknowns [kp_ij (j in nbr), kI_ij]"""
    blocks, rhs = [], []
    for r in range(mm):
        d = D[r]; zj = d['z'][nbr]; Cc = np.hstack([d['lam'] * zj, zj])[None, :]; Cr, br = realify(Cc, np.array([d['b'][i]])); blocks.append(Cr); rhs.append(br)
    return np.vstack(blocks), np.concatenate(rhs)


def lsq(C, b, rtol=1e-12):
    """min-norm LS via SVD with column-normalised rank decision; returns x, relative residual, rank"""
    sc = np.linalg.norm(C, axis=0); sc[sc == 0] = 1; Cn = C / sc
    U, s, Vt = np.linalg.svd(Cn, full_matrices=False); r = int(np.sum(s > rtol * s[0])) if s.size and s[0] > 0 else 0
    x = (Vt[:r].T @ ((U[:, :r].T @ b) / s[:r])) / sc
    res = np.linalg.norm(b - C @ x) / np.linalg.norm(b); return x, res, r


def free_nhop(D, mm, dist, n):
    """returns per-site eps, rank list, and full Kp, KI (absolute min-norm k_i) as 10x10"""
    nb = neighbourhoods(dist, n); eps = np.zeros(10); ranks = []; Kp = np.zeros((10, 10)); KI = np.zeros((10, 10))
    for i in range(10):
        C, b = rowsys(D, mm, i, nb[i]); x, res, rk = lsq(C, b); eps[i] = res; ranks.append(rk); k = len(nb[i]); Kp[i, nb[i]] = x[:k]; KI[i, nb[i]] = x[k:]
    return eps, ranks, Kp, KI


def poly_systems(D, mm, Lap, n, node_varying):
    """(b) shared / (c) node-varying graph polynomials. returns per-site relative residual (global/local LS)."""
    Lp = [np.linalg.matrix_power(Lap, l).astype(float) for l in range(n + 1)]
    def row_cols(i):   # matrix of size (2m) x 2(n+1): columns a_0..a_n, b_0..b_n, for site i
        cols = []
        for kind in (0, 1):
            for l in range(n + 1):
                col = []
                for r in range(mm):
                    d = D[r]; zl = Lp[l][i] * d['z']; val = (d['lam'] if kind == 0 else 1.0) * zl.sum(); col += [val.real, val.imag]
                col = np.array(col); cols.append(col)
        return np.array(cols).T
    bi = lambda i: np.concatenate([[D[r]['b'][i].real, D[r]['b'][i].imag] for r in range(mm)])
    # fix: rhs ordering must match row_cols ordering (re,im per target) -> consistent
    if node_varying:
        eps = np.array([lsq(row_cols(i), bi(i))[1] for i in range(10)])
    else:
        C = np.vstack([row_cols(i) for i in range(10)]); b = np.concatenate([bi(i) for i in range(10)]); x, gres, rk = lsq(C, b)
        eps = np.array([np.linalg.norm(bi(i) - row_cols(i) @ x) / np.linalg.norm(bi(i)) for i in range(10)])
    return eps


if __name__ == '__main__':
    m, tg = targets(5); D = modal_data(m, tg)
    tinfo = [dict(m=r + 1, root_re=d['root'].real, root_im=d['root'].imag, f_hz=d['root'].imag / 2 / np.pi, pll_participation=d['part'], lam_target=str(d['lam']), norm_q_over_null=float(la.norm(d['q'])),
                  norm_y=float(la.norm(d['y']))) for r, d in enumerate(D)]
    pd.DataFrame(tinfo).to_csv(CAMP / 'derived' / 'E1_targets.csv', index=False); print(pd.DataFrame(tinfo).to_string(), flush=True)
    rows, nmin, arch = [], [], []
    for gname, br in (('frozen', False), ('bridge', True)):
        g = comm_graph(br); dist = g['dist']; diam = g['diam']
        for mm in range(1, 6):
            ok = None
            for n in range(0, diam + 1):
                eps, ranks, _, _ = free_nhop(D, mm, dist, n)
                rows.append(dict(graph=gname, m=mm, n=n, eps=eps.max(), argmax_site=30 + int(eps.argmax()), min_rank=min(ranks), n_eq_per_site=2 * mm, n_unk_min=2 * min(len(x) for x in neighbourhoods(dist, n)),
                                 **{f'eps_site{30 + i}': eps[i] for i in range(10)}))
                if ok is None and eps.max() <= TOL: ok = n
                arch.append(dict(graph=gname, m=mm, n=n, arch='a_free', eps=eps.max()))
                for nm, nv in (('b_shared_poly', False), ('c_nodevarying_poly', True)):
                    e = poly_systems(D, mm, g['L'], n, nv); arch.append(dict(graph=gname, m=mm, n=n, arch=nm, eps=e.max()))
            nmin.append(dict(graph=gname, m=mm, n_min=(np.nan if ok is None else ok), eps_at_max_n=rows[-1]['eps'], largest_finite_n=diam))
    pd.DataFrame(rows).to_csv(CAMP / 'derived' / 'E1_residuals.csv', index=False); pd.DataFrame(nmin).to_csv(CAMP / 'derived' / 'E1_nmin.csv', index=False)
    pd.DataFrame(arch).to_csv(CAMP / 'derived' / 'E2_architectures.csv', index=False)
    print(pd.DataFrame(nmin).to_string()); print(pd.DataFrame(rows)[['graph', 'm', 'n', 'eps', 'min_rank', 'n_unk_min', 'n_eq_per_site']].to_string())
    print(pd.DataFrame(arch).pivot_table(index=['graph', 'm', 'n'], columns='arch', values='eps').to_string())
