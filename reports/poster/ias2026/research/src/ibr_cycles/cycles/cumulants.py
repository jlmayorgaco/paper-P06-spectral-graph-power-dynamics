from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..actions.green import ActionGreen

ComplexMatrix = NDArray[np.complex128]


def activation_potential(
    green: ActionGreen, theta: NDArray[np.float64], s: complex
) -> complex:
    """Phi(theta, s) = log det(I + Theta K(s)) with one scalar per action.

    Theta scales the rows of K belonging to each action, so the derivatives at
    theta = 0 are infinitesimal-amplitude quantities. They are NOT the same
    object as the finite-amplitude Moebius externality of mobius.py, and the
    two must never be plotted on the same axis without saying so.
    """

    k = green.k(s)
    scale = np.zeros(k.shape[0], dtype=np.complex128)
    for value, block in zip(theta, green.block_slices, strict=True):
        scale[block] = value
    identity = np.eye(k.shape[0], dtype=np.complex128)
    sign, magnitude = np.linalg.slogdet(identity + scale[:, None] * k)
    return complex(np.log(sign) + magnitude)


def second_cumulant(green: ActionGreen, a: int, b: int, s: complex) -> complex:
    """Analytic value of the mixed second derivative: -trace(K_ab K_ba)."""

    k = green.k(s)
    return -complex(np.trace(green.block(k, a, b) @ green.block(k, b, a)))


def third_cumulant(green: ActionGreen, a: int, b: int, c: int, s: complex) -> complex:
    """Analytic third cumulant: trace(K_ab K_bc K_ca) + trace(K_ac K_cb K_ba)."""

    k = green.k(s)

    def blk(i: int, j: int) -> ComplexMatrix:
        return green.block(k, i, j)

    forward = np.trace(blk(a, b) @ blk(b, c) @ blk(c, a))
    backward = np.trace(blk(a, c) @ blk(c, b) @ blk(b, a))
    return complex(forward + backward)


def finite_difference_cumulant(
    green: ActionGreen, axes: tuple[int, ...], s: complex, step: float
) -> complex:
    """Central mixed finite difference of Phi at theta = 0 along ``axes``."""

    n = len(green.factors)
    order = len(axes)
    total = 0.0 + 0.0j
    for mask in range(2**order):
        signs = [1 if (mask >> bit) & 1 else -1 for bit in range(order)]
        theta = np.zeros(n)
        for axis, sign in zip(axes, signs, strict=True):
            theta[axis] += sign * step
        weight = float(np.prod(signs))
        total += weight * activation_potential(green, theta, s)
    return total / (2.0 * step) ** order


@dataclass(frozen=True)
class CumulantConvergence:
    """Finite-difference convergence study against the analytic cumulant."""

    axes: tuple[int, ...]
    s: complex
    analytic: complex
    steps: tuple[float, ...]
    errors: tuple[float, ...]

    @property
    def best_error(self) -> float:
        return min(self.errors)

    @property
    def best_step(self) -> float:
        return self.steps[int(np.argmin(self.errors))]


def convergence_study(
    green: ActionGreen,
    axes: tuple[int, ...],
    s: complex,
    steps: tuple[float, ...] = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 5e-4, 2e-4),
) -> CumulantConvergence:
    """Compare the analytic cumulant with central differences over a step ladder."""

    if len(axes) == 2:
        analytic = second_cumulant(green, axes[0], axes[1], s)
    elif len(axes) == 3:
        analytic = third_cumulant(green, axes[0], axes[1], axes[2], s)
    else:
        raise ValueError("only second and third cumulants have closed forms here")
    scale = max(abs(analytic), 1e-300)
    errors = tuple(
        float(abs(finite_difference_cumulant(green, axes, s, step) - analytic) / scale)
        for step in steps
    )
    return CumulantConvergence(
        axes=axes, s=s, analytic=analytic, steps=steps, errors=errors
    )
