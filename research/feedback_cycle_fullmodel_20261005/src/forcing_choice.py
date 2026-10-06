"""Frozen rule (PREREGISTRATION): forcing bus = event bus in {8,16,29} with the largest linear input residue norm onto the limiting 44 ms root."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
J = CAMP / 'raw' / 'jac'; m = Model(P0, np.full(10, .044)); cat, _ = catalog(m, re_floor=-3, fmax=20)
lam, _ = m.refine(max([z for z in cat if z.imag > .3], key=lambda z: z.real)); r, l, _ = null_vectors(m, lam); den = np.vdot(l, m.derivative(lam) @ r)
x0 = pd.read_csv(J / 'x0_rho0p875.csv').x0.to_numpy(); rows = []
for bus in (8, 16, 29):
    b = pd.read_csv(J / f'inp_bus{bus}_rho0p875.csv').dx.to_numpy().astype(complex); dv = pd.read_csv(J / f'inpv_bus{bus}_rho0p875.csv').dv.to_numpy()
    for k, (th, _, _) in enumerate(pll_index()):
        ur, ui = dv[2 * (29 + k)], dv[2 * (29 + k) + 1]; e_in = -np.sin(x0[th]) * ur + np.cos(x0[th]) * ui
        b -= m.B[:, k] * e_in * (1 - np.exp(-lam * m.tau[k]))
    rows.append(dict(bus=bus, residue_norm=abs(np.vdot(l, b)) / abs(den), lam=str(lam), f_star_hz=lam.imag / 2 / np.pi))
T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE11a_forcing_bus_residue.csv', index=False); print(T.to_string()); print('chosen bus', int(T.loc[T.residue_norm.idxmax(), 'bus']), 'f* =', lam.imag / 2 / np.pi)
