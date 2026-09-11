"""BC01: a port representation that is regular at s = 0 (T1a after the quotient).

With the DAE blocks (f_x block-diagonal by device, f_z, g_x, g_z), the port
operator on the bus voltages and its identity are

    T(s) = g_z + g_x (sI - f_x)^-1 f_z,
    det P(s) = det(sI - f_x) det T(s) = det(g_z) det(sI - A),
    A = f_x - f_z g_z^-1 g_x.

Every portfolio of a governor-free, zero-damping model has det(sI - A) with a
structural double zero at s = 0: the rotation R (A R = 0) and its Jordan
partner w (A w = omega_B R). det T(s) inherits it. In floating point the double
zero splits by O(sqrt(Jacobian error)), so neither det T(0) nor the ratio
det T_S / det T_0 carries any information at the origin. That is why the
original port test could not certify an aperiodic crossing.

After the quotient, two exact relocations are applied in the DEVICE block, so
that the port Schur complement is still taken onto the network:

1. f_x -> f_x - beta U U^T with U = R/||R||. By T1a the rotation moves to -beta.
2. x0 = w + a R with a = omega_B/beta - U^T w/||R|| is an exact null vector of
   the relocated A. Then f_x -> f_x - beta2 x0 x0^T/||x0||^2, and by Brauer's
   theorem the neutral mode moves to -beta2. Nothing else moves.

The relocated port operator T##(s) is regular at 0 unless a physical eigenvalue
sits there:

    det T##(0) = det(g_z) det(-A##) / det(-f_x##).

A physical real eigenvalue crossing 0 flips sign(det(-A##)). The flip is
visible at the port unless the device factor det(-f_x##) flips too, which
would mean a device-internal mode. That factor is reported as the ledger.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..models.ieee39_devices import OMEGA_B


@dataclass(frozen=True)
class RelocatedPort:
    t0_sign: float  # sign of det T##(0)
    t0_logabs: float
    device_sign: float  # sign of det(-f_x##): the h ledger at s = 0
    network_sign: float  # sign of det(g_z)
    a_sign: float  # sign of det(-A##) (quotient-level truth)
    identity_residual: float  # |log|det T##(0)| - log|det g_z det(-A##)/det(-f_x##)||
    min_abs_eig_relocated: float


def _logdet(m: np.ndarray) -> tuple[float, float]:
    sign, logabs = np.linalg.slogdet(m)
    return float(sign), float(logabs)


def relocated_port(jac, r_x, w, beta: float = 1.0, beta2: float = 1.0) -> RelocatedPort:
    fx, fz, gx, gz = jac.fx, jac.fz, jac.gx, jac.gz
    r_norm = float(np.linalg.norm(r_x))
    u = r_x / r_norm
    fx1 = fx - beta * np.outer(u, u)
    alpha = OMEGA_B / beta - float(u @ w) / r_norm
    x0 = w + alpha * r_x
    fx2 = fx1 - beta2 * np.outer(x0, x0) / float(x0 @ x0)
    t0 = gz - gx @ np.linalg.solve(fx2, fz)  # T##(0) = g_z + g_x (0 - f_x##)^-1 f_z
    a2 = fx2 - fz @ np.linalg.solve(gz, gx)
    st, lt = _logdet(t0)
    sd, ld = _logdet(-fx2)
    sn, ln = _logdet(gz)
    sa, la = _logdet(-a2)
    return RelocatedPort(
        t0_sign=st,
        t0_logabs=lt,
        device_sign=sd,
        network_sign=sn,
        a_sign=sa,
        identity_residual=abs(lt - (ln + la - ld)),
        min_abs_eig_relocated=float(np.abs(np.linalg.eigvals(a2)).min()),
    )
