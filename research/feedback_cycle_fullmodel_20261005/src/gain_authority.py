"""Phase H: exact simple-root gain sensitivities d lambda/d(log Kp_i), d lambda/d(log Ki_i) from the full descriptor Delta(s) (exact exponentials)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
def sens(m, lam):
    """returns (dKp[10], dKi[10]) = d lambda / d(log gain) for the exact full characteristic; l^H Delta_p r / l^H Delta_s r."""
    r, l, _ = null_vectors(m, lam); e = np.exp(-lam * m.tau)
    den = np.vdot(l, (m.derivative(lam)) @ r)
    Ce = m.C @ r
    dKp = np.array([np.vdot(l, M['Bp'][:, i]) * e[i] * Ce[i] / den * m.p[10 + i] for i in range(10)])   # -(-l^H Bp e C r)/den * Kp
    dKi = np.array([np.vdot(l, M['Bi'][:, i]) * e[i] * Ce[i] / den * m.p[20 + i] for i in range(10)])
    return dKp, dKi
if __name__ == '__main__':
    m = Model(P0, np.full(10, .044)); cat, _ = catalog(m, re_floor=-3, fmax=20)
    crit = [z for z in cat if z.real > -SIGMA and z.imag > .3]; rows = []
    # finite-difference check of the sensitivity formula (first root, site 33)
    z0 = crit[0]; dKp, dKi = sens(m, z0); i = 3; eps = 1e-5
    pp = P0.copy(); pp[10 + i] *= np.exp(eps); pm = P0.copy(); pm[10 + i] *= np.exp(-eps)
    zp, _ = Model(pp, m.tau).refine(z0); zm, _ = Model(pm, m.tau).refine(z0); fd = (zp - zm) / (2 * eps)
    pp = P0.copy(); pp[20 + i] *= np.exp(eps); pm = P0.copy(); pm[20 + i] *= np.exp(-eps)
    zp2, _ = Model(pp, m.tau).refine(z0); zm2, _ = Model(pm, m.tau).refine(z0); fd2 = (zp2 - zm2) / (2 * eps)
    print('FD check Kp:', abs(fd - dKp[i]) / abs(fd), ' Ki:', abs(fd2 - dKi[i]) / abs(fd2), flush=True)
    for z in crit:
        z, _ = m.refine(z); dKp, dKi = sens(m, z)
        for i in range(10):
            rows.append(dict(root_f_hz=z.imag / 2 / np.pi, root_re=z.real, site=30 + i, dRe_dlogKp=dKp[i].real, dRe_dlogKi=dKi[i].real, dIm_dlogKp=dKp[i].imag, dIm_dlogKi=dKi[i].imag,
                             logKp_room_down=np.log(P0[10 + i] / KPLO), logKp_room_up=np.log(KPHI / P0[10 + i]), logKi_room_down=np.log(P0[20 + i] / KILO), logKi_room_up=np.log(KIHI / P0[20 + i])))
    T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE08_gain_authority.csv', index=False)
    lead = T[T.root_re == T.root_re.max()]
    print('leading root Re=%.4f f=%.3f Hz' % (lead.root_re.iloc[0], lead.root_f_hz.iloc[0]))
    print(lead[['site', 'dRe_dlogKp', 'dRe_dlogKi']].round(3).to_string(index=False))
    # P4 metric: best Kp-only vs best Ki-only |dRe/dlogK| (max over sites) for each critical root
    for f_, g in T.groupby('root_f_hz'):
        print('f=%.3f Hz Re=%+.3f | max|dRe/dlogKp|=%.3f  max|dRe/dlogKi|=%.3f  ratio=%.3f' % (f_, g.root_re.iloc[0], g.dRe_dlogKp.abs().max(), g.dRe_dlogKi.abs().max(), g.dRe_dlogKp.abs().max() / g.dRe_dlogKi.abs().max()))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True); roots = sorted(T.root_f_hz.unique())
    for ax, col, ttl in zip(axs, ['dRe_dlogKp', 'dRe_dlogKi'], ['d Re(lambda) / d log Kp_i [1/s]', 'd Re(lambda) / d log Ki_i [1/s]']):
        A = np.array([T[T.root_f_hz == f_].sort_values('site')[col].to_numpy() for f_ in roots]); v = np.abs(A).max()
        im = ax.imshow(A, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto'); ax.set_xticks(range(10)); ax.set_xticklabels(range(30, 40)); ax.set_yticks(range(len(roots))); ax.set_yticklabels([f'{f_:.2f} Hz' for f_ in roots]); ax.set_title(ttl, fontsize=9); ax.set_xlabel('PLL site (bus)')
        fig.colorbar(im, ax=ax)
    fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'FIG07_gain_authority_heatmap.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG07_gain_authority_heatmap.pdf')
