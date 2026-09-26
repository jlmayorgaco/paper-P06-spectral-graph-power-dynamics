from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .dae import LinearSystem, RealMatrix, RealVector, SemiExplicitDAE

_EPS = float(np.finfo(np.float64).eps)


@dataclass(frozen=True)
class DaeJacobians:
    """The four blocks of the linearized semi-explicit DAE."""

    fx: RealMatrix
    fz: RealMatrix
    gx: RealMatrix
    gz: RealMatrix

    @property
    def gz_condition(self) -> float:
        return float(np.linalg.cond(self.gz))


def _adaptive_step(value: float, scale: float) -> float:
    """Central-difference step of order eps^(1/3) times the variable scale."""

    return float(np.cbrt(_EPS) * max(abs(value), scale))


def central_difference_jacobians(
    model: SemiExplicitDAE,
    x: RealVector,
    z: RealVector,
    theta: dict[str, float],
    *,
    scale_x: float = 1.0,
    scale_z: float = 1.0,
) -> DaeJacobians:
    """Second-order accurate Jacobians of f and g at the given point."""

    x = np.asarray(x, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    n_x, n_z = x.size, z.size
    fx = np.zeros((n_x, n_x))
    gx = np.zeros((n_z, n_x))
    for j in range(n_x):
        h = _adaptive_step(float(x[j]), scale_x)
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        fx[:, j] = (model.f(xp, z, theta) - model.f(xm, z, theta)) / (2.0 * h)
        gx[:, j] = (model.g(xp, z, theta) - model.g(xm, z, theta)) / (2.0 * h)
    fz = np.zeros((n_x, n_z))
    gz = np.zeros((n_z, n_z))
    for j in range(n_z):
        h = _adaptive_step(float(z[j]), scale_z)
        zp, zm = z.copy(), z.copy()
        zp[j] += h
        zm[j] -= h
        fz[:, j] = (model.f(x, zp, theta) - model.f(x, zm, theta)) / (2.0 * h)
        gz[:, j] = (model.g(x, zp, theta) - model.g(x, zm, theta)) / (2.0 * h)
    return DaeJacobians(fx=fx, fz=fz, gx=gx, gz=gz)


def reduce_index_one(
    jac: DaeJacobians,
    labels: tuple[str, ...],
    *,
    condition_limit: float = 1e12,
) -> LinearSystem:
    """Return the reduced dynamic matrix fx - fz gz^-1 gx for nonsingular gz."""

    condition = jac.gz_condition
    if not np.isfinite(condition) or condition > condition_limit:
        raise ValueError(
            f"gz is numerically singular (cond={condition:.3e}); the model is "
            "not index-1 at this operating point"
        )
    a_red = jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)
    return LinearSystem(A=np.asarray(a_red, dtype=np.float64), labels=labels)


def jacobian_agreement(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
    """Relative Frobenius disagreement between two Jacobian estimates."""

    denominator = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1e-300)
    return float(np.linalg.norm(left - right) / denominator)
