"""Phase C: SG->GFL continuation (uniform rho 0->0.875, tau=40 ms, exact exponential delay). Branches tracked by Newton continuation + MAC of the
(right) null vectors, never by frequency alone. (1) PLL-family branches are tracked DOWN from rho=0.875; (2) all-SG branches are tracked UP from rho=0.125,
seeded by SG-subspace MAC against the native all-SG (rho=0, 114-state) eigenvectors."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
J = CAMP / 'raw' / 'jac'; g = groups_for()
def mac(a, b): return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))
rows = []
def record(branch, kind, r, z, res, v, m, mc):
    pp = participation(m, z, g)
    rows.append(dict(branch=branch, kind=kind, rho=r, real=z.real, imag=z.imag, f_hz=z.imag / 2 / np.pi, zeta=-z.real / abs(z), PLL_part=pp['PLL'], GFL_part=pp['GFL'], SG_part=pp['SG'], MAC_prev=mc, residual=res))
# ---- (1) PLL family down
m875 = Model(design(.875), TAU); cat, _ = catalog(m875, re_floor=-30, fmax=20)
fam = [z for z in cat if z.imag > 5 and participation(m875, z, g)['PLL'] >= .5]
print('PLL-family roots at rho=.875:', len(fam))
grid = np.round(np.arange(.875, .1249, -.0125), 4)
for bi, z0 in enumerate(fam):
    z = z0; vprev = None
    for r in grid:
        m = Model(design(r), TAU)
        try: z, res = m.refine(z)
        except Exception as e: print('branch', bi, 'lost at rho', r, e); break
        v, l, _ = null_vectors(m, z); mc = 1.0 if vprev is None else mac(vprev, v); vprev = v
        record(f'PLL{bi}', 'pll_family_down', r, z, res, v, m, mc)
# ---- (2) all-SG branches up
Fx = pd.read_csv(J / 'Fx_rho0p0.csv').to_numpy(); sm0 = pd.read_csv(J / 'statemap_rho0p0.csv'); sm1 = pd.read_csv(J / 'statemap_rho0p125.csv')
w, vl, vr = la.eig(Fx, right=True, left=True)
key0 = {(int(r.site), int(r.local_pos)): int(r['index']) - 1 for _, r in sm0[sm0.kind == 'SG'].iterrows()}
key1 = {(int(r.site), int(r.local_pos)): int(r['index']) - 1 for _, r in sm1[sm1.kind == 'SG'].iterrows()}
common = sorted(set(key0) & set(key1)); i0 = [key0[k] for k in common]; i1 = [key1[k] for k in common]
m125 = Model(design(.125), TAU); cat125, _ = catalog(m125, re_floor=-5, fmax=15)
vec125 = [null_vectors(m125, z)[0] for z in cat125]
osc = [k for k in np.argsort(-w.real) if w[k].imag > 0.1]; used = set(); unmatched = []
for bi, k in enumerate(osc):
    r0 = vr[:, k][i0]; best = max(range(len(cat125)), key=lambda j: mac(r0, vec125[j][i1]))
    mc0 = mac(r0, vec125[best][i1]); z = cat125[best]; vprev = None
    if mc0 < 0.7 or best in used: unmatched.append((bi, w[k], mc0)); continue
    used.add(best)
    rows.append(dict(branch=f'SG{bi}', kind='allSG_native', rho=0.0, real=w[k].real, imag=w[k].imag, f_hz=w[k].imag / 2 / np.pi, zeta=-w[k].real / abs(w[k]), PLL_part=np.nan, GFL_part=np.nan, SG_part=1.0, MAC_prev=np.nan, residual=np.nan))
    for r in np.round(np.arange(.125, .8751, .0125), 4):
        m = Model(design(r), TAU)
        try: z, res = m.refine(z)
        except Exception as e: print('SG branch', bi, 'lost at', r); break
        v, l, _ = null_vectors(m, z); mc = mc0 if vprev is None else mac(vprev, v); vprev = v
        record(f'SG{bi}', 'allSG_up', r, z, res, v, m, mc)
print('all-SG modes not matched at rho=0.125 (MAC<0.7 or duplicate):', [(b, round(z.imag/6.2832,4), round(mc,3)) for b, z, mc in unmatched])
T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE_C01_mode_continuation.csv', index=False)
fin = T[(T.rho == .875) & (T.kind != 'allSG_native')]
print(fin[['branch', 'kind', 'real', 'f_hz', 'zeta', 'PLL_part', 'MAC_prev']].round(4).to_string())
print('min MAC per branch:'); print(T.groupby('branch').MAC_prev.min().round(3).sort_values().head(12).to_string())
fig, ax = plt.subplots(figsize=(7.4, 4.6))
for b, gg in T.groupby('branch'):
    if gg.kind.iloc[0] == 'pll_family_down': ax.plot(gg.rho, gg.f_hz, color='#B3363A', lw=1.5)
    else: ax.plot(gg.rho, gg.f_hz, color='#176494', lw=1, alpha=.8)
ax.set_ylim(0, 8); ax.set_xlabel('uniform SG->GFL share rho'); ax.set_ylabel('frequency [Hz]'); ax.grid(alpha=.2)
ax.plot([], [], color='#B3363A', label='PLL-family branches (tracked from rho=0.875)'); ax.plot([], [], color='#176494', label='all-SG branches (tracked up from rho=0)'); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'FIG_C01_mode_migration.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG_C01_mode_migration.pdf')
