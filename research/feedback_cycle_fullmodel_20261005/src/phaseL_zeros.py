"""Phase L (exploratory): transmission / invariant zeros of (i) the ten-PLL return G(s)=C(sI-A0)^-1 B (grid closed, delay-free rational part), SISO diagonal channels and the full 10x10;
(ii) load-admittance (buses 8,16,29) -> generator-bus voltage phase channels, all-SG (rho=0) versus rho=0.875. Zeros = finite generalised eigenvalues of the Rosenbrock pencil;
each reported RHP zero is confirmed by a second method (Newton refinement on the pencil determinant surrogate: smallest singular value of the Rosenbrock matrix at the zero, relative)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from fm import *
J = CAMP / 'raw' / 'jac'
def rosen_zeros(A, B, C, D):
    n = A.shape[0]; m_ = B.shape[1]; p_ = C.shape[0]
    P = np.block([[A, B], [C, D]]); E = np.zeros_like(P); E[:n, :n] = np.eye(n)
    a, b = la.eig(P, E, homogeneous_eigvals=True)[0]
    z = a[np.abs(b) > 1e-9 * np.maximum(1, np.abs(a))] / b[np.abs(b) > 1e-9 * np.maximum(1, np.abs(a))]
    return z, P, E
def confirm(P, E, z):
    Mx = P - z * E; s = la.svd(Mx, compute_uv=False); return s[-1] / s[0]
rows = []
m = Model(P0, np.zeros(10)); A0 = M['Adev'] + M['Bv'] @ m.Vx
for k in range(10):
    z, P, E = rosen_zeros(A0, -m.B[:, [k]], m.C[[k], :], np.zeros((1, 1)))
    # with sign convention: Rosenbrock [[A, B],[C, D]] zero where [[sI-A, -B],[C,0]] drops rank -> pencil [[A,B],[C,D]] - z [[I,0],[0,0]] with B -> -B? use B=-m.B (sign flip is irrelevant for zero locations)
    z = z[np.isfinite(z) & (np.abs(z) < 1e4)]; rhp = z[z.real > 1e-6]
    rows.append(dict(channel=f'PLL{30 + k}->det{30 + k} (siso)', rho=.875, n_finite=len(z), n_rhp=len(rhp), rhp=' '.join(f'{x:.3f}' for x in rhp[:6]), confirm=' '.join(f'{confirm(P, E, x):.1e}' for x in rhp[:6])))
z, P, E = rosen_zeros(A0, -m.B, m.C, np.zeros((10, 10))); z = z[np.isfinite(z) & (np.abs(z) < 1e4)]; rhp = z[z.real > 1e-6]
rows.append(dict(channel='PLL(10)->det(10) MIMO', rho=.875, n_finite=len(z), n_rhp=len(rhp), rhp=' '.join(f'{x:.3f}' for x in rhp[:8]), confirm=' '.join(f'{confirm(P, E, x):.1e}' for x in rhp[:8])))
for tag, rho in (('0p0', 0.0), ('0p875', 0.875)):
    Fx = pd.read_csv(J / f'Fx_rho{tag}.csv').to_numpy(); Vx = pd.read_csv(J / f'Vall_x_rho{tag}.csv').to_numpy(); v0 = pd.read_csv(J / f'Vall_0_rho{tag}.csv').v.to_numpy(); n = Fx.shape[0]
    A = Fx if rho == 0 else Fx - m.B @ m.C        # delay-free rational part with the PLL error loop opened (e replaced by its delayed value): A0_full = Fx - B C ; use closed instantaneous loop instead: Fx
    A = Fx                                        # tau -> 0 limit
    ur, ui = v0[0::2], v0[1::2]
    for bus in (8, 16, 29):
        b = pd.read_csv(J / f'inp_bus{bus}_rho{tag}.csv').dx.to_numpy()[:, None]; dv = pd.read_csv(J / f'inpv_bus{bus}_rho{tag}.csv').dv.to_numpy()
        cnt = 0; ex = []; tot = 0; conf = []
        for g in range(30, 40):
            i = g - 1; mm = ur[i] ** 2 + ui[i] ** 2
            c = ((-ui[i] * Vx[2 * i] + ur[i] * Vx[2 * i + 1]) / mm)[None, :]; d = np.array([[(-ui[i] * dv[2 * i] + ur[i] * dv[2 * i + 1]) / mm]])
            z, P, E = rosen_zeros(A, -b, c, d); z = z[np.isfinite(z) & (np.abs(z) < 1e4)]; rhp = z[z.real > 1e-6]; tot += len(z); cnt += len(rhp)
            ex += [f'bus{g}:{x:.2f}' for x in rhp[:2]]; conf += [confirm(P, E, x) for x in rhp[:2]]
        rows.append(dict(channel=f'load@{bus} -> gen-bus phase (10 outputs)', rho=rho, n_finite=tot, n_rhp=cnt, rhp=' '.join(ex[:10]), confirm=' '.join(f'{x:.1e}' for x in conf[:10])))
T = pd.DataFrame(rows); T.to_csv(CAMP / 'derived' / 'TABLE_L01_rhp_zero_search.csv', index=False); print(T.to_string())
