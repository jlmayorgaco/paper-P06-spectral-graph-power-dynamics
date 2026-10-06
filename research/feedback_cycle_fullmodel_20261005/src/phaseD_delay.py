"""Phase D: exact-delay continuation tau 40->50 ms (uniform), baseline gains, rho=0.875. Tracks every root of the 44 ms catalogue with Re>-3 and
every PLL-family branch, by Newton continuation + MAC. Also N_unstable / N_margin by the banded argument-principle counter (floating point, NOT a certificate)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
from t3_full_spectrum_count import count
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
g = groups_for(); P = design(.875)
def mac(a, b): return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))
m44 = Model(P, np.full(10, .044)); cat44, _ = catalog(m44, re_floor=-6, fmax=20)
print('roots at 44 ms with Re>-0.05 (non-gauge):', int(np.sum(cat44.real > -SIGMA)), '| Re>0.01:', int(np.sum(cat44.real > .01)), flush=True)
# start tracking at tau=40 from the catalogue of the 44 ms roots (go back to 40 first), then forward
starts40 = []; taus = np.round(np.arange(40, 50.01, .1), 3); rows = []
for bi, z44 in enumerate(cat44):
    if z44.real < -3 or z44.imag < 0.3: continue
    z = z44; branch = {}
    # back to 40 ms
    for t in np.round(np.arange(44, 39.99, -.05), 3):
        m = Model(P, np.full(10, t / 1000))
        try: z, res = m.refine(z)
        except Exception: z = None; break
        branch[t] = z
    if z is None: continue
    zz = branch[40.0]; vprev = None
    if any(abs(zz - q) < 1e-6 for q in starts40): print('DUPLICATE branch at 40 ms (branch jump) for 44 ms root', z44, flush=True); continue
    starts40.append(zz)
    for t in taus:
        m = Model(P, np.full(10, t / 1000))
        try: zz, res = m.refine(zz)
        except Exception: break
        v, l, sv = null_vectors(m, zz); mc = 1.0 if vprev is None else mac(vprev, v); vprev = v
        pp = participation(m, zz, g)
        Dh = m.eye * 0 + (zz * np.eye(204) - M['Adev']); cnd = float(np.linalg.cond(Dh))
        rows.append(dict(branch=f'R{bi}', tau_ms=t, real=zz.real, imag=zz.imag, f_hz=zz.imag / 2 / np.pi, PLL_part=pp['PLL'], MAC_prev=mc, hidden_block_cond=cnd, residual=res,
                         Kp_lo=KPLO, Kp_hi=KPHI, Ki_lo=KILO, Ki_hi=KIHI))
T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE_D01_delay_continuation.csv', index=False)
cross = []
for b, gg in T.groupby('branch'):
    gg = gg.sort_values('tau_ms'); r = gg.real.to_numpy(); t = gg.tau_ms.to_numpy()
    def x(th):
        i = np.where((r[:-1] < th) & (r[1:] >= th))[0]
        return float(t[i[0]] + (th - r[i[0]]) / (r[i[0] + 1] - r[i[0]]) * (t[i[0] + 1] - t[i[0]])) if len(i) else np.nan
    cross.append(dict(branch=b, f_hz_at_40=gg.f_hz.iloc[0], f_hz_at_44=float(gg[gg.tau_ms == 44].f_hz.iloc[0]) if (gg.tau_ms == 44).any() else np.nan,
                      PLL_part=gg.PLL_part.mean(), tau_cross_margin_ms=x(-SIGMA), tau_cross_0_ms=x(0.0), re_at_44=float(gg[gg.tau_ms == 44].real.iloc[0]) if (gg.tau_ms == 44).any() else np.nan, min_MAC=gg.MAC_prev.min()))
C = pd.DataFrame(cross).sort_values('tau_cross_margin_ms'); C.to_csv(CAMP / 'derived' / 'TABLE_D01b_delay_crossings.csv', index=False); print(C.round(4).to_string())
cnt = []
for t in (40, 41, 42, 43, 44, 46, 48, 50):
    mm = Model(P, np.full(10, t / 1000)); nu, e1 = count(mm, .01); ng, e2 = count(mm, -SIGMA)
    cnt.append(dict(tau_ms=t, N_unstable=nu, N_margin_excl_gauge=ng - 1, max_noninteger_winding=max(e1, e2))); print(cnt[-1], flush=True)
pd.DataFrame(cnt).to_csv(CAMP / 'derived' / 'TABLE_D01c_root_counts.csv', index=False)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for b, gg in T.groupby('branch'): ax[0].plot(gg.tau_ms, gg.real, lw=1.3, label=f'{gg.f_hz.iloc[0]:.2f} Hz'); ax[1].plot(gg.real, gg.imag / 2 / np.pi, lw=1)
ax[0].axhline(0, color='k', lw=.8); ax[0].axhline(-SIGMA, color='#B3363A', ls='--', lw=1); ax[0].set_xlabel('uniform PLL delay [ms]'); ax[0].set_ylabel('Re(lambda) [1/s]'); ax[0].set_ylim(-3, 3); ax[0].grid(alpha=.2); ax[0].legend(fontsize=6, frameon=False, ncol=2)
ax[1].axvline(0, color='k', lw=.8); ax[1].set_xlabel('Re(lambda) [1/s]'); ax[1].set_ylabel('frequency [Hz]'); ax[1].set_xlim(-3, 3); ax[1].grid(alpha=.2)
fig.tight_layout(); fig.savefig(CAMP / 'figures' / 'FIG_D01_delay_root_locus.png', dpi=200); fig.savefig(CAMP / 'figures' / 'FIG_D01_delay_root_locus.pdf')
