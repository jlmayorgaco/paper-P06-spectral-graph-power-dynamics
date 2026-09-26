from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from ..dynamics.dae import LinearSystem

STATE_LABELS = ("delta1", "delta2", "omega1", "omega2", "z")
RETAINED = (0, 1, 2, 3)
HIDDEN = (4,)


@dataclass(frozen=True)
class ToyParameters:
    """Two-angle network with one hidden first-order controller state.

        delta_dot = omega
        omega_dot = -L delta - D omega + k [1, -1]^T z
        z_dot     = -a z + omega1 - omega2

    The controller state z is the only hidden coordinate, which makes the
    angle-level self-energy available in closed form and turns every algebraic
    identity of the package into a machine-precision regression test.
    """

    b01: float = 0.85
    b02: float = 3.20
    b12: float = 0.73
    d: float = 0.587
    a: float = 2.40
    k: float = 0.47

    @property
    def laplacian(self) -> NDArray[np.float64]:
        return np.array(
            [[self.b01 + self.b12, -self.b12], [-self.b12, self.b02 + self.b12]],
            dtype=np.float64,
        )


ACTIONS: dict[str, dict[str, float]] = {
    "A": {"b01": 3.80},
    "B": {"b02": 1.85},
    "C": {"k": 1.10},
}


def apply_actions(base: ToyParameters, members: tuple[str, ...]) -> ToyParameters:
    """Additive action amplitudes on the declared parameter coordinates."""

    deltas: dict[str, float] = {}
    for name in members:
        if name not in ACTIONS:
            raise KeyError(f"unknown action {name!r}")
        for field_name, amount in ACTIONS[name].items():
            deltas[field_name] = deltas.get(field_name, 0.0) + amount
    resolved = {key: getattr(base, key) + value for key, value in deltas.items()}
    return replace(base, **resolved)


def state_matrix(parameters: ToyParameters) -> NDArray[np.float64]:
    """Closed-form reduced dynamic matrix; the model has no algebraic variables."""

    a = np.zeros((5, 5), dtype=np.float64)
    a[0, 2] = 1.0
    a[1, 3] = 1.0
    a[2:4, 0:2] = -parameters.laplacian
    a[2, 2] = -parameters.d
    a[3, 3] = -parameters.d
    a[2, 4] = parameters.k
    a[3, 4] = -parameters.k
    a[4, 2] = 1.0
    a[4, 3] = -1.0
    a[4, 4] = -parameters.a
    return a


def linear_system(parameters: ToyParameters) -> LinearSystem:
    return LinearSystem(A=state_matrix(parameters), labels=STATE_LABELS)


def delta_a(base: ToyParameters, name: str) -> NDArray[np.float64]:
    """State-matrix update of one action, evaluated exactly."""

    return state_matrix(apply_actions(base, (name,))) - state_matrix(base)


def angle_operator(parameters: ToyParameters, s: complex) -> NDArray[np.complex128]:
    """Angle-level operator with the exact controller self-energy.

        T_delta(s) = (s^2 + D s) I + L - [k s / (s + a)] b b^T,  b = [1, -1]^T

    EXACT IDENTITY: det(sI - A) = (s + a) det(T_delta(s)).
    """

    b = np.array([1.0, -1.0])
    dressing = (parameters.k * s / (s + parameters.a)) * np.outer(b, b)
    return (
        (s * s + parameters.d * s) * np.eye(2, dtype=np.complex128)
        + parameters.laplacian
        - dressing
    )
