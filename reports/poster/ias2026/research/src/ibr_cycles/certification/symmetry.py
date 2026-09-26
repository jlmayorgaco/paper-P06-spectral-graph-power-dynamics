"""T1 / T1a: quotient of a declared physical symmetry and its relocation.

Theory: spectral_portfolio_theory_v1/TEORIA_Y_DEMOSTRACIONES.md, section 4.

Generic part (any matrix)
    quotient(A, R)          A_q = Z^T A Z for an orthonormal complement Z of im R,
                            with the block residuals that make (5)-(6) checkable
    relocate(A, R, beta)    A - beta U U^T (T1a)
    deflate_eigenvector     the same construction for a known exact eigenvector
                            (used for the neutral-frequency mode, which is a
                            declared physical mode, NOT a symmetry)

Model part (the IEEE-39 / Kundur / IEEE-68 phasor DAE of this repository)
    rotation_generator      derived from the equations: every device angle state
                            (machine rotor angle, converter PLL angle) shifts by
                            one, every bus voltage rotates, v -> j v. All other
                            states are written in a device-local rotating frame
                            and are invariant. Constant-power and constant-
                            impedance loads and the network are rotation-covariant.
    frequency_partner       the uniform frequency shift. With zero rotor damping,
                            constant mechanical power and a frequency-independent
                            network, x(t) = x0 + e w + e omega_B t R solves the
                            linearized DAE, so A w = omega_B R: a Jordan partner
                            of the symmetry. Its quotient image is an EXACT zero
                            eigenvector of A_q, i.e. a physical marginal mode that
                            must be reported, never deleted by size.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import null_space

from ..models.ieee39_devices import OMEGA_B

Matrix = NDArray[np.float64]


@dataclass(frozen=True)
class Quotient:
    """A_q and the evidence for the block-triangular identity (5)."""

    a_q: Matrix
    u: Matrix
    z: Matrix
    symmetry_residual: float  # ||A U|| / ||A||  (must be ~0: A R = 0)
    coupling_block: Matrix  # U^T A Z (the B_g block, may be nonzero)

    @property
    def dimension(self) -> int:
        return int(self.a_q.shape[0])


def orthonormal(r: Matrix) -> Matrix:
    r = np.atleast_2d(np.asarray(r, dtype=np.float64))
    if r.shape[0] == 1 and r.shape[1] > 1:
        r = r.T
    q, _ = np.linalg.qr(r)
    return q[:, : r.shape[1]]


def quotient(a: Matrix, r: Matrix) -> Quotient:
    """T1: [U Z]^T A [U Z] = [[0, B_g], [0, A_q]] when A R = 0."""

    a = np.asarray(a, dtype=np.float64)
    u = orthonormal(r)
    z = null_space(u.T)
    norm = max(float(np.linalg.norm(a, 2)), 1e-300)
    return Quotient(
        a_q=z.T @ a @ z,
        u=u,
        z=z,
        symmetry_residual=float(np.linalg.norm(a @ u, 2) / norm),
        coupling_block=u.T @ a @ z,
    )


def relocate(a: Matrix, r: Matrix, beta: float) -> Matrix:
    """T1a: shift only the declared symmetry subspace from 0 to -beta."""

    u = orthonormal(r)
    return np.asarray(a, dtype=np.float64) - beta * (u @ u.T)


def deflate_eigenvector(a_q: Matrix, v: Matrix) -> Quotient:
    """Quotient by an exact (right) eigenvector of eigenvalue zero.

    Identical algebra to T1: if A_q v = 0 then det(sI - A_q) = s det(sI - A_qq).
    The removed mode stays in the ledger of the caller.
    """

    return quotient(a_q, v)


# --------------------------------------------------------------------- model --


def _base(label: str) -> str:
    return label.rsplit("_", 1)[0]


def rotation_generator(dae, z: NDArray[np.float64]) -> tuple[Matrix, Matrix]:
    """(R_x, R_z): d/dphi of the global rotation, derived from the device frames.

    Rotor angles ``delta`` and PLL angles ``theta_pll`` shift by one; every other
    device state is written in its own rotating frame. Bus voltages rotate:
    d(v e^{j phi})/d phi = j v, i.e. (vx, vy) -> (-vy, vx).
    """

    r_x = np.zeros(dae.n_x)
    for k, label in enumerate(dae.labels):
        if _base(label) in ("delta", "theta_pll"):
            r_x[k] = 1.0
    z = np.asarray(z, dtype=np.float64)
    r_z = np.empty_like(z)
    r_z[0::2] = -z[1::2]
    r_z[1::2] = z[0::2]
    return r_x, r_z


@dataclass(frozen=True)
class FrequencyPartner:
    w: Matrix | None
    reason: str  # "DERIVED" or why none exists


def frequency_partner(dae) -> FrequencyPartner:
    """The uniform-frequency-shift direction w with A w = omega_B R_x.

    Derived device by device:
      two-axis machine (IEEE-39, Kundur): omega +1 (absolute speed, pu)
      68-bus machine: slip +1; speed-input PSS washout state +K (it tracks K*slip)
      GFL converter: PLL integrator x_pll +omega_B (steady drift with v_q = 0);
                     synthetic-inertia filter y_vi +1 (it tracks the PLL deviation)
    Exists only if every machine has zero rotor damping; otherwise returns None.
    """

    w = np.zeros(dae.n_x)
    for slot in dae.slots:
        device = slot.device
        p = getattr(device, "parameters", None)
        damping = getattr(p, "d", 0.0) if p is not None else 0.0
        if slot.kind == "sg" and damping != 0.0:
            return FrequencyPartner(
                None, f"rotor damping D={damping:g} at bus {slot.bus}"
            )
        for k, label in enumerate(device.labels):
            base = _base(label)
            index = slot.start + k
            if slot.kind == "sg" and base == "omega":
                w[index] = 1.0
            elif slot.kind == "sg" and base == "pss_w" and hasattr(p, "pss") and p.pss:
                w[index] = dict(p.pss)["k"]
            elif slot.kind == "gfl" and base == "x_pll":
                w[index] = OMEGA_B
            elif slot.kind == "gfl" and base == "y_vi":
                w[index] = 1.0
    return FrequencyPartner(w, "DERIVED")


@dataclass(frozen=True)
class DaeIdentities:
    """Residuals of the symmetry and Jordan identities on the DAE blocks."""

    f_rot: float  # ||f_x R_x + f_z R_z|| / scale
    g_rot: float  # ||g_x R_x + g_z R_z|| / scale
    a_rot: float  # ||A R_x|| / (||A|| ||R_x||)
    a_jordan: float | None  # ||A w - omega_B R_x|| / (||A|| ||w||)
    g_jordan: float | None  # ||g_x w|| / scale (frequency shift leaves injections)


def dae_identities(jac, a: Matrix, r_x, r_z, w) -> DaeIdentities:
    sx = max(np.linalg.norm(jac.fx, 2), np.linalg.norm(jac.fz, 2))
    sg = max(np.linalg.norm(jac.gx, 2), np.linalg.norm(jac.gz, 2))
    na = max(float(np.linalg.norm(a, 2)), 1e-300)
    rn = float(np.linalg.norm(np.concatenate([r_x, r_z])))
    out = DaeIdentities(
        f_rot=float(np.linalg.norm(jac.fx @ r_x + jac.fz @ r_z) / (sx * rn)),
        g_rot=float(np.linalg.norm(jac.gx @ r_x + jac.gz @ r_z) / (sg * rn)),
        a_rot=float(np.linalg.norm(a @ r_x) / (na * np.linalg.norm(r_x))),
        a_jordan=None,
        g_jordan=None,
    )
    if w is not None:
        nw = float(np.linalg.norm(w))
        out = DaeIdentities(
            f_rot=out.f_rot,
            g_rot=out.g_rot,
            a_rot=out.a_rot,
            a_jordan=float(np.linalg.norm(a @ w - OMEGA_B * r_x) / (na * nw)),
            g_jordan=float(np.linalg.norm(jac.gx @ w) / (sg * nw)),
        )
    return out
