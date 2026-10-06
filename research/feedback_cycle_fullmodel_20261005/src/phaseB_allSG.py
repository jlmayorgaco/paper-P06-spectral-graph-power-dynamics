"""Phase B: all-SG reference (rho=0, native 114-state support, no GFL/PLL states) modes + load-to-bus-frequency response."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
J = CAMP / 'raw' / 'jac'
Fx = pd.read_csv(J / 'Fx_rho0p0.csv').to_numpy(); sm = pd.read_csv(J / 'statemap_rho0p0.csv'); n = Fx.shape[0]
w, vl, vr = la.eig(Fx, left=True, right=True)
rows = []
for k in np.argsort(-w.real):
    z = w[k]
    if z.imag < 0.1: continue
    p = np.abs(vl[:, k].conj() * vr[:, k]); p /= p.sum(); site = {int(s): float(p[(sm.site == s).to_numpy()].sum()) for s in range(30, 40)}
    top = sorted(site, key=site.get, reverse=True)[:3]
    ang = np.angle(vr[:, k][sm[sm.local_pos == 1].index.to_numpy()])  # first SG state per site (rotor angle) phase pattern
    rows.append(dict(real=z.real, imag=z.imag, f_hz=z.imag / 2 / np.pi, zeta=-z.real / abs(z), PLL_participation='N/A', SG_participation=float(p.sum()),
                     top_sites=','.join(map(str, top)), top_site_weights=','.join(f'{site[t]:.3f}' for t in top)))
T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE03_allSG_modes.csv', index=False)
print(T.round(4).to_string(), flush=True)
print('max real part non-gauge:', sorted(w.real)[-3:])
# --- GFL baseline (rho=.875, tau=40ms) PLL-family modes for comparison
m = Model(P0, TAU); cat, _ = catalog(m, re_floor=-30, fmax=20); g = groups_for()
G = []
for z in cat:
    if z.imag < 0.5: continue
    pp = participation(m, z, g); G.append(dict(real=z.real, imag=z.imag, f_hz=z.imag / 2 / np.pi, zeta=-z.real / abs(z), PLL_participation=pp['PLL'], GFL_participation=pp['GFL'], SG_participation=pp['SG']))
GT = pd.DataFrame(G); GT.to_csv(CAMP / 'derived' / 'TABLE04a_baseline_GFL_modes_40ms.csv', index=False)
print(GT.round(4).head(30).to_string())
# --- response: |s * theta_bus| per MW of load admittance at buses 8,16,29; max over gen buses
def theta_rows(tag):
    Vx = pd.read_csv(J / f'Vall_x_rho{tag}.csv').to_numpy(); v0 = pd.read_csv(J / f'Vall_0_rho{tag}.csv').v.to_numpy()
    ur, ui = v0[0::2], v0[1::2]; rows_ = {}
    for b in range(30, 40):
        i = b - 1; m2 = ur[i] ** 2 + ui[i] ** 2
        rows_[b] = (-ui[i] * Vx[2 * i] + ur[i] * Vx[2 * i + 1]) / m2, (-ui[i], ur[i], m2, 2 * i)
    return rows_
f = np.geomspace(0.05, 20, 600); out = []
for tag, rho in (('0p0', 0.0), ('0p875', 0.875)):
    cth = theta_rows(tag); x0 = pd.read_csv(J / f'x0_rho{tag}.csv').x0.to_numpy()
    for bus in (8, 16, 29):
        b = pd.read_csv(J / f'inp_bus{bus}_rho{tag}.csv').dx.to_numpy(); dv = pd.read_csv(J / f'inpv_bus{bus}_rho{tag}.csv').dv.to_numpy()
        if rho == 0: A0 = Fx
        else: mm = Model(design(rho), TAU)
        pk = []
        for fi in f:
            s = 2j * np.pi * fi
            if rho == 0: X = la.solve(s * np.eye(n) - Fx, b)
            else:
                bs = b.astype(complex)
                # delayed PLL error: replace instantaneous e by e(t-tau) -> input correction  -B D_u (1 - e^{-s tau}) u  (D_u u from voltage feedthrough at PLL bus)
                for k, (th, _, _) in enumerate(pll_index()):
                    bus_k = 30 + k; ur_, ui_ = dv[2 * (bus_k - 1)], dv[2 * (bus_k - 1) + 1]
                    e_in = -np.sin(x0[th]) * ur_ + np.cos(x0[th]) * ui_
                    bs -= mm.B[:, k] * e_in * (1 - np.exp(-s * mm.tau[k]))
                X = la.solve(mm.delta(s), bs)
            th_ = [abs(s * (cth[bb][0] @ X + (cth[bb][1][0] * dv[cth[bb][1][3]] + cth[bb][1][1] * dv[cth[bb][1][3] + 1]) / cth[bb][1][2])) for bb in range(30, 40)]
            pk.append(max(th_) / (2 * np.pi))   # Hz per MW (0.1 MW basis already per MW)
        out.append(pd.DataFrame(dict(f_hz=f, bus=bus, rho=rho, max_gen_bus_freq_dev_Hz_per_MW=pk)))
FR = pd.concat(out); FR.to_csv(CAMP / 'derived' / 'TABLE03b_frf_load_to_busfreq.csv', index=False)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axs[0]; ax.scatter(w.real, w.imag / 2 / np.pi, s=14, c='#176494', label='all-SG (rho=0)')
ax.scatter(cat.real, cat.imag / 2 / np.pi, s=14, c='#C99A20', label='rho=0.875, 40 ms (exact delay)')
ax.set_xlim(-12, 0.3); ax.set_ylim(0, 8); ax.axvline(-0.05, color='#B3363A', ls='--', lw=1); ax.set_xlabel('Re(lambda) [1/s]'); ax.set_ylabel('frequency [Hz]'); ax.legend(frameon=False, fontsize=8); ax.grid(alpha=.2)
ax = axs[1]
for (rho, bus), g_ in FR.groupby(['rho', 'bus']):
    ax.loglog(g_.f_hz, g_.max_gen_bus_freq_dev_Hz_per_MW, ls='-' if rho == 0 else '--', label=f'rho={rho}, bus {bus}')
ax.set_xlabel('forcing frequency [Hz]'); ax.set_ylabel('max gen-bus |dF| [Hz per MW]'); ax.legend(frameon=False, fontsize=7); ax.grid(alpha=.2, which='both')
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'FIG_B01_allSG_spectrum.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG_B01_allSG_spectrum.pdf')
FR['pk'] = FR.max_gen_bus_freq_dev_Hz_per_MW
print(FR.groupby(['rho', 'bus']).apply(lambda g_: g_.loc[g_.pk.idxmax(), ['f_hz', 'pk']]).to_string())
