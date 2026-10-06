"""Replaced generation vs stability margin: rightmost non-gauge root (exact delay) vs uniform inverter share rho, tau = 40 and 44 ms, baseline gains."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
rows = []
for rho in [.5, .6, .7, .75, .8, .85, .875, .9, .925, .95]:
    for tms in (40, 44):
        m = Model(design(rho), np.full(10, tms / 1000)); cat, _ = catalog(m, re_floor=-3, fmax=20); z = max(cat, key=lambda q: q.real)
        rows.append(dict(rho=rho, tau_ms=tms, rightmost_re=z.real, f_hz=z.imag / 2 / np.pi, n_beyond_margin_pairs=int(sum(q.real > -SIGMA and q.imag > .3 for q in cat))))
        print(rows[-1], flush=True); pd.DataFrame(rows).to_csv(CAMP / 'derived' / 'TABLE13_share_margin.csv', index=False)
