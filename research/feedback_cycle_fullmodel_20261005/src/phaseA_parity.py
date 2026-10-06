"""Phase A: parity of the Python exact-delay reader with fresh Julia Jacobians; baseline roots and counts."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
from t3_full_spectrum_count import count
rows = []
for tag, r in [('0p125', .125), ('0p25', .25), ('0p375', .375), ('0p5', .5), ('0p625', .625), ('0p75', .75), ('0p875', .875)]:
    Fx = pd.read_csv(CAMP / 'raw' / 'jac' / f'Fx_rho{tag}.csv').to_numpy()
    m = Model(design(r), np.zeros(10)); err = la.norm(Fx - m.A) / la.norm(Fx)
    ev = la.eigvals(Fx); evp = la.eigvals(m.A)
    rows.append(dict(test='Fx_julia_vs_python_rel_fro', rho=r, value=err, tol=1e-8, passed=bool(err < 1e-8),
                     max_eig_diff=float(np.max(np.abs(np.sort_complex(ev) - np.sort_complex(evp))))))
    print(rows[-1], flush=True)
# baseline catalogued roots (historical TABLE_02) refined with the exact characteristic
cat = pd.read_csv(ROOT / 'experiments/graph_gsp_codesign_20261003/TABLE_02_BASELINE_ROOTS.csv')
m = Model(P0, TAU); d = []
for z in cat.real.to_numpy() + 1j * cat.imag.to_numpy():
    r, e = m.refine(z); d.append(abs(r - z))
rows.append(dict(test='baseline_11_roots_refine_shift', rho=.875, value=max(d), tol=1e-8, passed=bool(max(d) < 1e-8)))
# own fresh catalogue vs historical: every historical root recovered?
cr, ce = catalog(m, re_floor=-5, fmax=15)
miss = [z for z in cat.real.to_numpy() + 1j * cat.imag.to_numpy() if min(abs(cr - z)) > 1e-6]
rows.append(dict(test='fresh_catalog_recovers_historical_roots', rho=.875, value=len(miss), tol=0, passed=len(miss) == 0))
for tms in (40, 44):
    mm = Model(P0, np.full(10, tms / 1000)); nu, e1 = count(mm, 0.01); ng, e2 = count(mm, -SIGMA)
    rows.append(dict(test=f'full_count_tau{tms}ms', rho=.875, value=nu, N_unstable=nu, N_margin_excl_gauge=ng - 1, max_noninteger_winding=max(e1, e2), passed=True))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(CAMP / 'derived' / 'TABLE02_baseline_parity.csv', index=False)
print(pd.DataFrame(rows).to_string())
