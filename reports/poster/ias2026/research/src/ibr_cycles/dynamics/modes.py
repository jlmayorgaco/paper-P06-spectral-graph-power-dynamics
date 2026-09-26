from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class Mode:
    """One eigenvalue with the diagnostics needed to trust it."""

    value: complex
    right: NDArray[np.complex128]
    left: NDArray[np.complex128]
    condition: float
    participation: NDArray[np.float64]

    @property
    def real(self) -> float:
        return float(self.value.real)

    @property
    def imag(self) -> float:
        return float(self.value.imag)

    @property
    def frequency_hz(self) -> float:
        return float(abs(self.value.imag) / (2.0 * np.pi))

    @property
    def damping(self) -> float:
        magnitude = abs(self.value)
        if magnitude == 0.0:
            return 1.0
        return float(-self.value.real / magnitude)


@dataclass(frozen=True)
class Spectrum:
    """All modes of a reduced dynamic matrix, sorted by decreasing real part."""

    modes: tuple[Mode, ...]

    @property
    def values(self) -> NDArray[np.complex128]:
        return np.array([mode.value for mode in self.modes], dtype=np.complex128)

    @property
    def critical(self) -> Mode:
        return self.modes[0]

    @property
    def spectral_abscissa(self) -> float:
        return self.critical.real

    @property
    def is_stable(self) -> bool:
        return self.spectral_abscissa < 0.0

    def rhp_count(self, *, boundary: float = 0.0) -> int:
        return int(np.sum(self.values.real > boundary))


def eigen_analysis(matrix: NDArray[np.float64] | ComplexMatrix) -> Spectrum:
    """Eigenvalues with left/right vectors, condition numbers and participations."""

    a = np.asarray(matrix, dtype=np.complex128)
    values, right = np.linalg.eig(a)
    left = np.linalg.inv(right).conj().T
    order = np.argsort(-values.real, kind="stable")
    modes: list[Mode] = []
    for index in order:
        r = right[:, index]
        left_vector = left[:, index]
        overlap = complex(left_vector.conj() @ r)
        condition = float(
            np.linalg.norm(left_vector) * np.linalg.norm(r) / max(abs(overlap), 1e-300)
        )
        weights = np.abs(left_vector.conj() * r)
        total = float(weights.sum())
        participation = weights / total if total > 0.0 else weights
        modes.append(
            Mode(
                value=complex(values[index]),
                right=r,
                left=left_vector,
                condition=condition,
                participation=np.asarray(participation, dtype=np.float64),
            )
        )
    return Spectrum(modes=tuple(modes))
