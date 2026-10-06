import pathlib, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
CAMP = pathlib.Path(__file__).resolve().parents[1]; DER = CAMP / 'derived'; FIG = CAMP / 'figures'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 18, 'axes.labelsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'legend.fontsize': 14, 'pdf.fonttype': 42})
SZ = (136 / 25.4, 84 / 25.4); COL = ['#003B2D', '#B3363A', '#176494', '#C99A20', '#6A3D9A']; FLOOR = 1e-16
def fl(x): return np.maximum(x, FLOOR)
def new():
    fig, ax = plt.subplots(figsize=SZ); ax.set_yscale('log'); ax.grid(alpha=.25, which='major'); return fig, ax
E1 = pd.read_csv(DER / 'E1_residuals.csv'); A = pd.read_csv(DER / 'E2_architectures.csv'); E5 = pd.read_csv(DER / 'E5_spillover.csv')
for g, suf in (('frozen', ''), ('bridge', '_bridge')):
    fig, ax = new()
    for m in range(1, 6):
        d = E1[(E1.graph == g) & (E1.m == m)]; ax.plot(d.n, fl(d.eps), 'o-', color=COL[m - 1], label=f'm={m}', lw=2, ms=6)
    ax.axhline(1e-8, color='k', ls=':', lw=1.5); ax.set_xlabel('hops n'); ax.set_ylabel(r'$\epsilon(n)$'); ax.set_ylim(5e-17, 1); ax.set_xticks(sorted(E1[E1.graph == g].n.unique()))
    ax.legend(ncol=3, frameon=False, loc='upper right', fontsize=12, columnspacing=.8, handlelength=1.2); fig.tight_layout(pad=.3); fig.savefig(FIG / f'E1_residual_vs_hops{suf}.pdf'); plt.close(fig)
    fig, ax = new(); m = 4; st = {'a_free': ('o-', COL[0], 'free n-hop'), 'b_shared_poly': ('s--', COL[1], 'shared poly'), 'c_nodevarying_poly': ('^-.', COL[2], 'node poly')}
    for k, (mk, c, lab) in st.items():
        d = A[(A.graph == g) & (A.m == m) & (A.arch == k)]; ax.plot(d.n, fl(d.eps), mk, color=c, label=lab, lw=2, ms=6)
    ax.axhline(1e-8, color='k', ls=':', lw=1.5); ax.set_xlabel('hops / polynomial degree n'); ax.set_ylabel(r'$\epsilon(n)$, m=4'); ax.set_ylim(5e-17, 1); ax.set_xticks(sorted(A[A.graph == g].n.unique()))
    ax.legend(frameon=False, loc='lower left', fontsize=12); fig.tight_layout(pad=.3); fig.savefig(FIG / f'E2_architectures{suf}.pdf'); plt.close(fig)
    fig, ax = plt.subplots(figsize=SZ); ax.grid(alpha=.25)
    for m in range(1, 6):
        for var, mk, fc in (('delta', 'o-', None), ('absolute', 's--', 'none')):
            d = E5[(E5.graph == g) & (E5.m == m) & (E5.variant == var)]
            if len(d): ax.plot(d.n, d.N_margin, mk, color=COL[m - 1], lw=1.5, ms=6, mfc=(fc or COL[m - 1]), label=f'm={m}' if var == 'delta' else None)
    ax.axhline(10, color='k', ls=':', lw=1.5); ax.axhline(0, color='k', lw=.8); ax.set_xlabel('hops n'); ax.set_ylabel('roots beyond -0.05'); ax.set_xticks(sorted(E1[E1.graph == g].n.unique()))
    ax.legend(ncol=3, frameon=False, loc='upper right', fontsize=11, columnspacing=.8, handlelength=1.2); ax.set_ylim(-1, 24); fig.tight_layout(pad=.3); fig.savefig(FIG / f'E5_spillover{suf}.pdf'); plt.close(fig)
