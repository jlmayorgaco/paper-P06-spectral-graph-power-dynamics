from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .self_energy import SelfEnergy, StatePartition


@dataclass(frozen=True)
class SchurResidual:
    """Evidence for the exact identity det(sI-A) = det(sI-Ahh) det(T_eff(s))."""

    points: tuple[complex, ...]
    relative_errors: tuple[float, ...]

    @property
    def max_relative_error(self) -> float:
        return max(self.relative_errors)


def schur_determinant_residual(
    matrix: NDArray[np.generic],
    partition: StatePartition,
    points: Sequence[complex],
) -> SchurResidual:
    """Verify the Schur determinant factorization pointwise.

    EXACT IDENTITY. Any deviation beyond roundoff is an implementation bug, not
    a modelling choice, so this is a hard gate rather than a diagnostic.
    """

    a = np.asarray(matrix, dtype=np.complex128)
    n = a.shape[0]
    sigma = SelfEnergy.from_matrix(a, partition)
    errors: list[float] = []
    for s in points:
        full = np.linalg.det(s * np.eye(n, dtype=np.complex128) - a)
        hidden = np.linalg.det(
            s * np.eye(sigma.n_hidden, dtype=np.complex128) - sigma.ahh
        )
        effective = np.linalg.det(sigma.effective_operator(s))
        scale = max(abs(full), abs(hidden * effective), 1e-300)
        errors.append(float(abs(full - hidden * effective) / scale))
    return SchurResidual(points=tuple(points), relative_errors=tuple(errors))
