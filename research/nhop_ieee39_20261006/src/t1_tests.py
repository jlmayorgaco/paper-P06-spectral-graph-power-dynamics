"""T1 tests: parity, counter validation, open-PLL operator + determinant identity, graph export."""
import sys, pathlib, json, time
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from nhop import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng = np.random.default_rng(20261006); R = {}
# --- parity
mods = []
errs = []
for tau in (.040, .044):
    ref = Model(P0, np.full(10, tau)); dm = diag_model(tau)
    for _ in range(5):
        s = complex(rng.uniform(-1, 1), rng.uniform(0, 80))
        a, b = ref.delta(s), dm.delta(s); errs.append(la.norm(a - b) / la.norm(a))
R['parity_max_rel_err'] = float(max(errs)); R['parity_pass'] = bool(max(errs) <= 1e-12); print('parity', R['parity_max_rel_err'], flush=True)
# --- root refinement vs fm.catalog
m44 = diag_model(.044); ref44 = Model(P0, np.full(10, .044)); cat, _ = catalog(ref44, re_floor=-3, fmax=20)
crit = [z for z in cat if z.real > -SIGMA and z.imag > 0.3]; d = []
for z in crit:
    z2, r = m44.refine(z); z1, _ = ref44.refine(z); d.append(abs(z2 - z1))
R['refine_vs_Model_max_diff'] = float(max(d)); print('refine diff', R['refine_vs_Model_max_diff'], flush=True)
# --- counter
t0 = time.time()
for ms in (40, 44):
    nu, nm, w = spectrum_counts(diag_model(ms * 1e-3)); R[f'count_{ms}ms'] = dict(N_unstable=nu, N_margin=nm, max_noninteger_winding=w, seconds=time.time() - t0)
    print(ms, R[f'count_{ms}ms'], flush=True)
# --- open PLL operator and determinant identity
op = OpenPLL(); R['open_pll_checks'] = op.chk; print(op.chk)
ids = []
for k in range(24):
    s = complex(rng.uniform(-0.5, 2), rng.uniform(0.5, 90)); tau = rng.uniform(.01, .06)
    Kp = rng.normal(size=(10, 10)) * 20; KI = rng.normal(size=(10, 10)) * 200
    dm = DModel(Kp, KI, tau)
    sl, ld = np.linalg.slogdet(dm.delta(s))
    sr, lr = np.linalg.slogdet(s * np.eye(op.nr) - op.Arr)
    E = np.diag(np.exp(-s * dm.tau)); Mx = s**2 * (1 + TF * s) * np.eye(10) - (s * Kp + KI) @ E @ op.G(s)
    s2, l2 = np.linalg.slogdet(Mx)
    # c = tf^-10
    ratio = (sl / (sr * s2)) * np.exp(ld - lr - l2) / TF**-10
    ids.append(abs(ratio - 1))
R['identity_constant'] = 'c = t_f^-10 = (600 pi)^10'; R['identity_max_rel_err'] = float(max(ids)); R['identity_pass'] = bool(max(ids) <= 1e-8)
print('identity', R['identity_max_rel_err'], flush=True)
# --- graph
for tag, br in (('frozen', False), ('bridge', True)):
    g = comm_graph(br); A = g['A']; edges = [(30 + i, 30 + j, g['R'][i, j]) for i in range(10) for j in range(i + 1, 10) if A[i, j]]
    pd.DataFrame(edges, columns=['bus_i', 'bus_j', 'r_ij']).assign(graph=tag).to_csv(CAMP / 'derived' / ('GRAPH_Gc.csv' if tag == 'frozen' else 'GRAPH_Gc_bridge.csv'), index=False)
    pd.DataFrame(g['dist']).to_csv(CAMP / 'derived' / f'GRAPH_Gc_{tag}_hopdist.csv', index=False)
    R[f'graph_{tag}'] = dict(n_edges=len(edges), components=g['info']['components'], labels=g['info']['labels'], largest_finite_diam=g['diam'], bridge=g['info']['bridge_added'],
                             edges=[(a, b) for a, b, _ in edges], degrees=A.sum(1).tolist())
    if tag == 'frozen': R['Y_block_convention'] = g['info']
# figure (layout: classical MDS of electrical distance)
g = comm_graph(True); Rm = g['R']; J = np.eye(10) - 1 / 10; Bm = -.5 * J @ (Rm**2) @ J; w, V = np.linalg.eigh(Bm); pos = V[:, -2:] * np.sqrt(np.maximum(w[-2:], 0))
fig, ax = plt.subplots(figsize=(5.2, 4.2)); gf = comm_graph(False)['A']
for i in range(10):
    for j in range(i + 1, 10):
        if g['A'][i, j]: ax.plot(*pos[[i, j]].T, color='#B3363A' if not gf[i, j] else '#003B2D', ls='--' if not gf[i, j] else '-', lw=2)
for i in range(10): ax.text(*pos[i], str(30 + i), ha='center', va='center', fontsize=10, bbox=dict(boxstyle='circle', fc='#FFF4D6', ec='#003B2D'))
ax.set_title('G_c (green: frozen rule; red dashed: AMENDMENT_01 bridge)', fontsize=8); ax.axis('off'); ax.set_aspect('equal'); fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'Gc.pdf')
json.dump(R, open(CAMP / 'derived' / 'T1_tests.json', 'w'), indent=1, default=str)
