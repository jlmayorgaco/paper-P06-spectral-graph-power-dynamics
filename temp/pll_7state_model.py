"""
PLL 7-state IBR model + second-order damping-margin correction.
Versioned reference implementation for cross-verification (Claude <-> ChatGPT).

Model: 3 nodes. N0 = synchronous generator (SG), N2 = grid-forming (GFM),
N1 = grid-following (GFL) inverter with a 2nd-order PLL (states: angle, freq,
integral z1).

States x = [theta0, theta1_pll, theta2, omega0, omega_pll1, omega2, z1]
- theta_i : bus angles (nodes 0,1,2)
- omega_i : node frequencies (0,1,2). For the GFL node, omega_pll1 is the PLL
  frequency state.
- z1      : PLL integral state (this is what breaks the clean 2nd-order/QEP form)

A_sys structure (7x7):
  rows 0,1,2 : d(theta_i)/dt = omega_i
  row 3 (SG)  : M0 * d(omega0)/dt = -L[0,:].theta - D0 omega0
  row 5 (GFM) : M2 * d(omega2)/dt = -L[2,:].theta - D2 omega2
  row 4 (GFL/PLL): d(omega_pll1)/dt = -Kp*(L[1,:].theta) + Ki*z1 - dp*omega_pll1
  row 6 (PLL integral): d(z1)/dt = -L[1,:].theta

NOTE: this is a *physically-motivated reduced* IBR model for method development,
NOT a substitute for ANDES REGCA1/REGCP1. The PLL adds a genuine state (z1),
so A_sys is non-normal and is NOT a clean 2nd-order QEP.
"""
import numpy as np


def build_L(a, b, c):
    """3-node weighted Laplacian with edge weights a=(0-1), b=(1-2), c=(0-2)."""
    return np.array([[a + c, -a, -c],
                     [-a, a + b, -b],
                     [-c, -b, b + c]], float)


def build_gfl_7state(L, M0, M2, Kp, Ki, dp, D0=0.4, D2=0.5):
    """Assemble the 7x7 state matrix of the SG+GFM+GFL(PLL) system."""
    A = np.zeros((7, 7))
    A[0, 3] = 1.0          # dtheta0 = omega0
    A[1, 4] = 1.0          # dtheta1 = omega_pll1
    A[2, 5] = 1.0          # dtheta2 = omega2
    A[3, 0:3] = -L[0, :] / M0
    A[3, 3] = -D0 / M0     # SG swing
    A[5, 0:3] = -L[2, :] / M2
    A[5, 5] = -D2 / M2     # GFM swing
    A[4, 0:3] = -Kp * L[1, :]
    A[4, 6] = Ki
    A[4, 4] = -dp          # PLL frequency dynamics
    A[6, 0:3] = -L[1, :]   # PLL integral state
    return A


def zeta_min_full(A):
    """Exact damping margin from the full eigenvalues of A_sys."""
    ev = np.linalg.eigvals(A)
    zs = [-l.real / abs(l) for l in ev
          if abs(l.imag) > 1e-6 and abs(l) > 1e-9]
    return min(zs) if zs else 1.0


def zeta_second_order(A, off_state=(4, 6)):
    """
    Second-order perturbation estimate of the damping margin for the
    first-order system A_sys.

    A0 = A with the PLL integral coupling switched off (off_state set to 0):
         this is the 'order-0' approximation (PLL as nodal damping, no memory).
    Delta A = A - A0 is the coupling that the PLL integral introduces.

    The critical pole of A0 is corrected to 1st and 2nd order using
    standard non-symmetric eigenvalue perturbation with left/right eigenvectors.

    Returns (zeta_order0, zeta_order1, zeta_order2).
    """
    A0 = A.copy()
    A0[off_state] = 0.0
    DA = A - A0

    ev0, V0 = np.linalg.eig(A0)
    W0 = np.linalg.inv(V0).conj().T   # rows are left eigenvectors

    # pick critical mode of A0 (smallest damping ratio among oscillatory modes)
    cand = [(-l.real / abs(l), i) for i, l in enumerate(ev0)
            if abs(l.imag) > 1e-6 and abs(l) > 1e-9]
    if not cand:
        return 1.0, 1.0, 1.0
    _, ic = min(cand, key=lambda t: t[0])
    lc, vc, wc = ev0[ic], V0[:, ic], W0[:, ic]

    # first-order correction
    l1 = (wc.conj() @ DA @ vc) / (wc.conj() @ vc)

    # second-order correction (sum over other modes)
    l2 = 0.0
    for l in range(len(ev0)):
        if l == ic:
            continue
        den = (wc.conj() @ vc) * (W0[:, l].conj() @ V0[:, l]) * (lc - ev0[l])
        if abs(den) > 1e-12:
            l2 += (wc.conj() @ DA @ V0[:, l]) * (W0[:, l].conj() @ DA @ vc) / den

    z = lambda s: -s.real / abs(s)
    return z(lc), z(lc + l1), z(lc + l1 + l2)


if __name__ == "__main__":
    # Canonical reproducible example
    A = build_gfl_7state(build_L(5, 4, 3), M0=10, M2=6, Kp=1.0, Ki=0.5, dp=2.0)
    zr = zeta_min_full(A)
    z0, z1, z2 = zeta_second_order(A)
    print(f"zeta_min (full 7-state) = {zr:.5f}")
    print(f"order-0  = {z0:.5f}  err {abs(z0-zr)/zr*100:.2f}%")
    print(f"order-1  = {z1:.5f}  err {abs(z1-zr)/zr*100:.2f}%")
    print(f"order-2  = {z2:.5f}  err {abs(z2-zr)/zr*100:.2f}%")

    # Monte Carlo robustness (matches the numbers reported to ChatGPT)
    rng = np.random.default_rng(11)
    e0, e1, e2 = [], [], []
    for _ in range(400):
        L = build_L(*rng.uniform(2, 8, 3))
        A = build_gfl_7state(L, rng.uniform(4, 12), rng.uniform(3, 8),
                             rng.uniform(0.5, 2), rng.uniform(0.2, 0.8),
                             rng.uniform(1, 3))
        if np.max(np.linalg.eigvals(A).real) > 1e-6:
            continue
        zr = zeta_min_full(A)
        if zr < 1e-4 or zr > 0.9:
            continue
        z0, z1, z2 = zeta_second_order(A)
        e0.append(abs(z0-zr)/zr); e1.append(abs(z1-zr)/zr); e2.append(abs(z2-zr)/zr)
    e0, e1, e2 = map(np.array, (e0, e1, e2))
    print(f"\nMonte Carlo ({len(e0)} stable cases):")
    print(f"order-0 : median {np.median(e0)*100:.2f}%  p95 {np.percentile(e0,95)*100:.1f}%")
    print(f"order-1 : median {np.median(e1)*100:.2f}%  p95 {np.percentile(e1,95)*100:.1f}%")
    print(f"order-2 : median {np.median(e2)*100:.2f}%  p95 {np.percentile(e2,95)*100:.1f}%")
