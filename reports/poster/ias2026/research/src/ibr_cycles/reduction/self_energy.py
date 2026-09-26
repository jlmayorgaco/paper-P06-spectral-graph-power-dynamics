from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class StatePartition:
    """Split of the reduced state vector into retained and hidden coordinates."""

    retained: tuple[int, ...]
    hidden: tuple[int, ...]

    @classmethod
    def from_retained(cls, n: int, retained: tuple[int, ...]) -> StatePartition:
        kept = tuple(sorted(retained))
        if len(set(kept)) != len(kept) or any(i < 0 or i >= n for i in kept):
            raise ValueError("retained indices must be unique and in range")
        hidden = tuple(i for i in range(n) if i not in set(kept))
        return cls(retained=kept, hidden=hidden)


@dataclass(frozen=True)
class SelfEnergy:
    """Frequency-dependent controller self-energy of the hidden states.

        Sigma(s) = Arh (sI - Ahh)^-1 Ahr

    This is the object that dresses propagation through the retained
    coordinates. It is never replaced by a constant nodal damping matrix: the
    whole mechanism under study lives in its frequency dependence.
    """

    arr: ComplexMatrix
    arh: ComplexMatrix
    ahr: ComplexMatrix
    ahh: ComplexMatrix

    @classmethod
    def from_matrix(
        cls, matrix: NDArray[np.generic], partition: StatePartition
    ) -> SelfEnergy:
        a = np.asarray(matrix, dtype=np.complex128)
        r = list(partition.retained)
        h = list(partition.hidden)
        return cls(
            arr=a[np.ix_(r, r)],
            arh=a[np.ix_(r, h)],
            ahr=a[np.ix_(h, r)],
            ahh=a[np.ix_(h, h)],
        )

    @property
    def n_retained(self) -> int:
        return int(self.arr.shape[0])

    @property
    def n_hidden(self) -> int:
        return int(self.ahh.shape[0])

    def __call__(self, s: complex) -> ComplexMatrix:
        if self.n_hidden == 0:
            return np.zeros((self.n_retained, self.n_retained), dtype=np.complex128)
        identity = np.eye(self.n_hidden, dtype=np.complex128)
        return self.arh @ np.linalg.solve(s * identity - self.ahh, self.ahr)

    def effective_operator(self, s: complex) -> ComplexMatrix:
        """T_eff(s) = sI - Arr - Sigma(s)."""

        identity = np.eye(self.n_retained, dtype=np.complex128)
        return s * identity - self.arr - self(s)

    def norm(self, s: complex) -> float:
        return float(np.linalg.norm(self(s), 2))
