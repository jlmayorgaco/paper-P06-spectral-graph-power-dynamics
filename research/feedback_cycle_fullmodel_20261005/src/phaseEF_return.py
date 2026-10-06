"""Phases E-G: feedback return G(s), locally dressed interaction Q(s) (AMENDMENT_01), parity, closure, cycles, SCCs, physical-graph comparison."""
import sys, pathlib, itertools
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng = np.random.default_rng(20261005); EYE = np.eye(204)
def parts(m, s):
    D = s * EYE - M['Adev'] - M['Bv'] @ m.Vx          # sI - A0
    G = m.C @ la.solve(D, m.B); e = np.exp(-s * m.tau)
    L = 1 - np.diag(G) * e
    Q = -(G - np.diag(np.diag(G))) * e[None, :] / L[:, None]
    return D, G, e, L, Q
# ---- parity: det Delta / (det(sI-A0) prod L) = det(I+Q)
m40 = Model(P0, TAU); m44 = Model(P0, np.full(10, .044)); rows = []
for k in range(24):
    m = m40 if k % 2 == 0 else m44; s = complex(rng.uniform(-0.3, 0.5), 2 * np.pi * rng.uniform(.5, 12))
    D, G, e, L, Q = parts(m, s)
    sl1, l1 = np.linalg.slogdet(m.delta(s)); sl2, l2 = np.linalg.slogdet(D); sl3, l3 = np.linalg.slogdet(np.eye(10) + Q)
    lhs = sl1 / (sl2 * np.prod(L / abs(L)) * sl3) * np.exp(l1 - l2 - np.sum(np.log(abs(L))) - l3)
    rows.append(dict(s=str(s), tau_ms=m.tau[0] * 1e3, parity_ratio_minus_1=abs(lhs - 1), passed=bool(abs(lhs - 1) < 1e-8)))
PAR = pd.DataFrame(rows); PAR.to_csv(CAMP / 'derived' / 'TABLE_E00_determinant_parity.csv', index=False)
print('determinant parity:', int(PAR.passed.sum()), '/', len(PAR), 'max |ratio-1| =', PAR.parity_ratio_minus_1.max(), flush=True)
assert PAR.passed.all(), 'sign convention failed parity: STOP'
# ---- roots to analyse
cat44, _ = catalog(m44, re_floor=-3, fmax=20); crit = [z for z in cat44 if z.real > -SIGMA and z.imag > 0.3]
targets = [('tau44', m44, z) for z in crit] + [('tau40', m40, complex(-0.8015607744782098, 30.850471693237324))]
edges, pairs, sccs, trip = [], [], [], []; ROOTQ = {}
for name, m, z in targets:
    z, res = m.refine(z); D, G, e, L, Q = parts(m, z); ev = la.eigvals(Q); k = np.argmin(abs(ev + 1)); sv = la.svd(np.eye(10) + Q, compute_uv=False)
    tag = f'{name}_{z.imag / 2 / np.pi:.3f}Hz'; ROOTQ[tag] = (z, Q, L)
    thr = .1 * np.max(abs(Q)); adj = (abs(Q.T) >= thr)       # edge i->j  if |Q_ji|>=thr  (information i feeds j)
    ncomp, lab = connected_components(adj, directed=True, connection='strong')
    for c in range(ncomp):
        mem = [30 + i for i in np.where(lab == c)[0]]; sccs.append(dict(root=tag, scc=c, size=len(mem), buses=' '.join(map(str, mem)), edge_threshold=thr))
    for i in range(10):
        for j in range(10):
            if i != j: edges.append(dict(root=tag, src=30 + j, dst=30 + i, abs_Q=abs(Q[i, j]), phase_deg=np.degrees(np.angle(Q[i, j]))))
    for i, j in itertools.combinations(range(10), 2):
        p = Q[i, j] * Q[j, i]; pairs.append(dict(root=tag, Re=z.real, f_hz=z.imag / 2 / np.pi, bus_i=30 + i, bus_j=30 + j, abs_p=abs(p), phase_deg=np.degrees(np.angle(p)), abs_1mp=abs(1 - p), abs_log_1mp=abs(np.log(abs(1 - p)))))
    for i, j, k2 in itertools.combinations(range(10), 3):
        for (a, b, c) in ((i, j, k2), (i, k2, j)):
            p = Q[a, b] * Q[b, c] * Q[c, a]; trip.append(dict(root=tag, buses=f'{30 + a}-{30 + b}-{30 + c}', abs_p=abs(p), abs_1mp=abs(1 - p), phase_deg=np.degrees(np.angle(p))))
    print(f'{tag}: Re={z.real:+.4f} | min|1+eig(Q)|={abs(ev[k] + 1):.2e} | min sigma(I+Q)={sv[-1]:.2e} | min|L_i|={abs(L).min():.3f} max|L_i|={abs(L).max():.3f} | rho(Q)={max(abs(ev)):.3f} | SCC sizes {sorted([s_["size"] for s_ in sccs if s_["root"] == tag])}', flush=True)
pd.DataFrame(edges).to_csv(CAMP / 'derived' / 'TABLE_F01_feedback_edges.csv', index=False)
PR = pd.DataFrame(pairs); PR.to_csv(CAMP / 'derived' / 'TABLE_F02_pair_cycles.csv', index=False)
pd.DataFrame(sccs).to_csv(CAMP / 'derived' / 'TABLE_F03_sccs.csv', index=False)
TR = pd.DataFrame(trip); TR.sort_values('abs_1mp').groupby('root').head(10).to_csv(CAMP / 'derived' / 'TABLE_F04_triple_cycles_top10.csv', index=False)
# ---- physical graph (Kron port coupling)
Y = pd.read_csv(ROOT / 'experiments/graph_gsp_codesign_20261003/model/Y.csv').to_numpy()
Ph = np.zeros((10, 10))
for i in range(10):
    for j in range(10):
        if i != j: Ph[i, j] = np.linalg.norm(Y[2 * i:2 * i + 2, 2 * j:2 * j + 2]) / np.sqrt(np.linalg.norm(Y[2 * i:2 * i + 2, 2 * i:2 * i + 2]) * np.linalg.norm(Y[2 * j:2 * j + 2, 2 * j:2 * j + 2]))
phys = pd.DataFrame([dict(bus_i=30 + i, bus_j=30 + j, phys_coupling=Ph[i, j]) for i, j in itertools.combinations(range(10), 2)]); phys['phys_rank'] = phys.phys_coupling.rank(ascending=False)
phys.to_csv(CAMP / 'derived' / 'TABLE_G01_physical_coupling.csv', index=False)
FROZEN = {(35, 36), (30, 37), (33, 34)}; summ = []
for tag, g in PR.groupby('root'):
    g = g.merge(phys, on=['bus_i', 'bus_j']); g['dyn_rank'] = g.abs_log_1mp.rank(ascending=False); g['dyn_rank_closure'] = g.abs_1mp.rank()
    top = g.sort_values('abs_log_1mp', ascending=False).head(5); rho_s = spearmanr(g.phys_coupling, g.abs_p).statistic
    fr = {f'{a}-{b}': (int(g[(g.bus_i == a) & (g.bus_j == b)].dyn_rank.iloc[0]), int(g[(g.bus_i == a) & (g.bus_j == b)].phys_rank.iloc[0])) for a, b in sorted(FROZEN)}
    summ.append(dict(root=tag, top5_pairs=' | '.join(f'{r.bus_i}-{r.bus_j} (|1-p|={r.abs_1mp:.3f})' for r in top.itertuples()), frozen_pair_ranks_dyn_phys=str(fr), spearman_physcoupling_vs_abs_p=rho_s,
                     physical_top5=' | '.join(f'{r.bus_i}-{r.bus_j}' for r in g.sort_values('phys_rank').head(5).itertuples())))
SM = pd.DataFrame(summ); SM.to_csv(CAMP / 'derived' / 'TABLE_G02_dynamic_vs_physical.csv', index=False); print(SM.to_string())
# ---- sigma_min(I+Q) on the imaginary axis (FIG04)
f = np.geomspace(.05, 20, 4000); out = {}
for name, m in (('tau40', m40), ('tau44', m44)):
    vals = []
    for fi in f:
        _, G, e, L, Q = parts(m, 2j * np.pi * fi); vals.append((la.svd(np.eye(10) + Q, compute_uv=False)[-1], abs(L).min()))
    out[name] = np.array(vals)
pd.DataFrame({'f_hz': f, 'smin_IQ_40': out['tau40'][:, 0], 'minL_40': out['tau40'][:, 1], 'smin_IQ_44': out['tau44'][:, 0], 'minL_44': out['tau44'][:, 1]}).to_csv(CAMP / 'derived' / 'TABLE_E01_sigma_min_scan.csv', index=False)
fig, ax = plt.subplots(figsize=(7.4, 4.2))
ax.loglog(f, out['tau40'][:, 0], color='#003B2D', label='sigma_min(I+Q), 40 ms'); ax.loglog(f, out['tau44'][:, 0], color='#B3363A', label='sigma_min(I+Q), 44 ms')
ax.loglog(f, out['tau44'][:, 1], color='#B3363A', ls=':', label='min_i |L_i|, 44 ms (local factors)'); ax.loglog(f, out['tau40'][:, 1], color='#003B2D', ls=':', label='min_i |L_i|, 40 ms')
ax.set_xlabel('frequency [Hz] on s = j omega'); ax.set_ylabel('singular value'); ax.legend(frameon=False, fontsize=8); ax.grid(alpha=.2, which='both'); fig.tight_layout()
fig.savefig(CAMP / 'figures' / 'FIG04_sigma_min_IplusQ.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG04_sigma_min_IplusQ.pdf')
# ---- core graph and cycle plot for the rightmost 44 ms root
lead = max((k for k in ROOTQ if k.startswith('tau44')), key=lambda k: ROOTQ[k][0].real); z, Q, L = ROOTQ[lead]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8)); ax = axs[0]; ang = np.linspace(0, 2 * np.pi, 11)[:-1] + np.pi / 2; pos = np.c_[np.cos(ang), np.sin(ang)]
mx = np.max(abs(Q))
for i in range(10):
    for j in range(10):
        if i != j and abs(Q[i, j]) >= .1 * mx: ax.annotate('', xy=pos[i] * .93, xytext=pos[j] * .93, arrowprops=dict(arrowstyle='-|>', lw=3 * abs(Q[i, j]) / mx, color='#B3363A', alpha=.55, shrinkA=14, shrinkB=14, connectionstyle='arc3,rad=.12'))
for i in range(10): ax.text(*pos[i], str(30 + i), ha='center', va='center', fontsize=12, bbox=dict(boxstyle='circle', fc='#FFF4D6', ec='#003B2D'))
ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'|Q_ij| >= 0.1 max, root {lead}', fontsize=9)
ax = axs[1]; g = PR[PR.root == lead]
ax.scatter(g.phase_deg, g.abs_p, c=np.where([(a, b) in FROZEN for a, b in zip(g.bus_i, g.bus_j)], '#C99A20', '#176494'), s=40)
for r in g.sort_values('abs_p', ascending=False).head(6).itertuples(): ax.annotate(f'{r.bus_i}-{r.bus_j}', (r.phase_deg, r.abs_p), fontsize=8)
ax.axhline(1, color='k', lw=.6); ax.set_xlabel('phase of p_ij = Q_ij Q_ji [deg]'); ax.set_ylabel('|p_ij|'); ax.set_yscale('log'); ax.grid(alpha=.2); ax.set_title('pair cycles (gold = frozen 35-36, 30-37, 33-34)', fontsize=9)
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'FIG05_feedback_core_graph.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG05_feedback_core_graph.pdf')
fig.savefig(CAMP / 'figures' / 'FIG06_pair_cycle_magnitude_phase.png', dpi=200)
