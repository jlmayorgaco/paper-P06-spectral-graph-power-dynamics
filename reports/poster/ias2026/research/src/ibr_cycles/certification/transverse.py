"""Transverse (relative) stability: quotient by the exact structural center subspace.

For the frozen phasor DAE (D = 0 on every machine, no governor, a
frequency-independent network), the linearization A has the exact invariant
subspace

    C = span{R_x, w},   A R_x = 0,   A w = omega_B R_x,

where R_x is the rotational (reference-angle) gauge symmetry and w the physical
common-frequency drift. In the basis (R_x, w) the restriction of A to C is the
Jordan block [[0, omega_B], [0, 0]]: a chain of length two, not two independent
zeros. The two directions have different physical status:

    R_x  gauge: the equations are invariant under a rigid rotation (no physics)
    w    physical marginal mode: no primary frequency restoration in the model

With U an orthonormal basis of C and Z an orthonormal basis of its complement,
invariance of C gives Z^T A U = 0 and

    [U Z]^T A [U Z] = [[U^T A U, U^T A Z], [0, A_perp]],   A_perp = Z^T A Z,
    det(sI - A) = s^2 det(sI - A_perp).

A_perp is the dynamics on the quotient space X / C, the TRANSVERSE dynamics.
Transverse stability := alpha(A_perp) < 0. No eigenvalue is removed by size.

When the partner does not exist (a machine with D != 0, or governors), C reduces
to span{R_x} and A_perp is the plain T1 quotient.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import null_space

from ..models.ieee39_devices import OMEGA_B
from .symmetry import orthonormal


@dataclass(frozen=True)
class CenterSubspace:
    dimension: int
    z_c: np.ndarray  # right basis [R_x, w] (or [R_x])
    w_c: np.ndarray | None  # left basis with W_c^T Z_c = I (None if not computed)
    jordan: np.ndarray  # A restricted to C in the basis z_c
    invariance_residual: float  # ||A Z_c - Z_c J|| / (||A|| ||Z_c||)
    rotation_residual: float  # ||A R_x|| / (||A|| ||R_x||)
    left_residual: float | None  # ||W_c^T A - J W_c^T|| / ||A||
    left_condition: float | None

    @property
    def is_jordan_chain(self) -> bool:
        return self.dimension == 2 and abs(self.jordan[0, 1]) > 0.5 * OMEGA_B


def center_subspace(a, r_x, w=None, *, left: bool = True) -> CenterSubspace:
    a = np.asarray(a, dtype=np.float64)
    norm = max(float(np.linalg.norm(a, 2)), 1e-300)
    cols = [r_x] if w is None else [r_x, w]
    z_c = np.column_stack(cols)
    j, *_ = np.linalg.lstsq(z_c, a @ z_c, rcond=None)
    inv = float(np.linalg.norm(a @ z_c - z_c @ j) / (norm * np.linalg.norm(z_c)))
    rot = float(np.linalg.norm(a @ r_x) / (norm * np.linalg.norm(r_x)))
    w_c = left_res = left_cond = None
    if left:
        # the left generalized null space of the zero block: null((A^k)^T), k = dim C
        power = np.linalg.matrix_power(a, z_c.shape[1])
        u_, s_, vt_ = np.linalg.svd(power.T)
        cand = vt_[-z_c.shape[1] :].T  # right singular vectors of A^k^T, smallest
        gram = cand.T @ z_c
        left_cond = float(np.linalg.cond(gram))
        w_c = cand @ np.linalg.inv(gram).T  # so that W_c^T Z_c = I
        left_res = float(np.linalg.norm(w_c.T @ a - j @ w_c.T) / norm)
    return CenterSubspace(
        dimension=z_c.shape[1],
        z_c=z_c,
        w_c=w_c,
        jordan=j,
        invariance_residual=inv,
        rotation_residual=rot,
        left_residual=left_res,
        left_condition=left_cond,
    )


@dataclass(frozen=True)
class Transverse:
    a_perp: np.ndarray
    z: np.ndarray
    u: np.ndarray
    coupling_residual: float  # ||Z^T A U|| / ||A||: must vanish (C invariant)


def transverse_operator(a, r_x, w=None) -> Transverse:
    a = np.asarray(a, dtype=np.float64)
    u = orthonormal(np.column_stack([r_x] if w is None else [r_x, w]))
    z = null_space(u.T)
    norm = max(float(np.linalg.norm(a, 2)), 1e-300)
    return Transverse(
        a_perp=z.T @ a @ z,
        z=z,
        u=u,
        coupling_residual=float(np.linalg.norm(z.T @ a @ u, 2) / norm),
    )


def transverse_spectrum(a, r_x, w=None) -> np.ndarray:
    return np.linalg.eigvals(transverse_operator(a, r_x, w).a_perp)
