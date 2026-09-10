from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .self_energy import SelfEnergy

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class ControllerPole:
    """One hidden-controller pole of Sigma(s) with its matrix residue."""

    pole: complex
    residue: ComplexMatrix
    condition: float

    @property
    def frequency_hz(self) -> float:
        return float(abs(self.pole.imag) / (2.0 * np.pi))

    @property
    def damping(self) -> float:
        magnitude = abs(self.pole)
        return 1.0 if magnitude == 0.0 else float(-self.pole.real / magnitude)

    @property
    def residue_norm(self) -> float:
        return float(np.linalg.norm(self.residue, 2))

    def contribution(self, s: complex) -> ComplexMatrix:
        return self.residue / (s - self.pole)

    def detuning(self, target: complex) -> float:
        """Distance from a network mode to this controller pole."""

        return float(abs(target - self.pole))


@dataclass(frozen=True)
class PoleExpansion:
    """Partial-fraction expansion Sigma(s) = sum_j R_j / (s - p_j).

    Valid only when Ahh is diagonalizable. A hidden block built from cascaded
    controller integrators is typically NOT: an integrator feeding another
    integrator is a Jordan block at the origin, and no simple-pole expansion
    exists. ``status`` reports that case as ``DEFECTIVE_HIDDEN_BLOCK`` instead of
    returning a plausible but meaningless residue table.
    """

    poles: tuple[ControllerPole, ...]
    max_reconstruction_error: float
    diagonalizable: bool

    @property
    def status(self) -> str:
        if not self.diagonalizable:
            return "DEFECTIVE_HIDDEN_BLOCK"
        if self.max_reconstruction_error > 1e-8:
            return "EXPANSION_UNRESOLVED"
        return "SUCCESS"

    @property
    def trustworthy(self) -> bool:
        return self.status == "SUCCESS"

    def dominant(self, s: complex) -> ControllerPole:
        return max(self.poles, key=lambda p: np.linalg.norm(p.contribution(s), 2))

    def evaluate(self, s: complex) -> ComplexMatrix:
        total = sum(pole.contribution(s) for pole in self.poles)
        return np.asarray(total, dtype=np.complex128)


def expand_self_energy(
    sigma: SelfEnergy,
    *,
    probe_points: tuple[complex, ...] = (0.5 + 1.0j, -0.3 + 4.0j, 1.7 - 2.5j),
    condition_limit: float = 1e10,
) -> PoleExpansion:
    """Pole/residue decomposition of the controller self-energy."""

    values, right = np.linalg.eig(sigma.ahh)
    condition = float(np.linalg.cond(right))
    diagonalizable = condition < condition_limit
    left = np.linalg.inv(right)
    poles: list[ControllerPole] = []
    for j, value in enumerate(values):
        residue = np.outer(sigma.arh @ right[:, j], left[j, :] @ sigma.ahr)
        poles.append(
            ControllerPole(
                pole=complex(value),
                residue=np.asarray(residue, dtype=np.complex128),
                condition=condition,
            )
        )
    expansion = PoleExpansion(
        poles=tuple(poles),
        max_reconstruction_error=np.inf,
        diagonalizable=diagonalizable,
    )
    errors = []
    for s in probe_points:
        exact = sigma(s)
        scale = max(float(np.linalg.norm(exact)), 1e-300)
        errors.append(float(np.linalg.norm(expansion.evaluate(s) - exact) / scale))
    return PoleExpansion(
        poles=tuple(poles),
        max_reconstruction_error=max(errors),
        diagonalizable=diagonalizable,
    )
