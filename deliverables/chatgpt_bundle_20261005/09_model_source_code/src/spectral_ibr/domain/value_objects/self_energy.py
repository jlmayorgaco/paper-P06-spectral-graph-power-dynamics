from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class SelfEnergy:
    """Self-energy Sigma(s) = D0 + C (sI - Acc)^-1 B."""

    d0: ComplexMatrix
    c: ComplexMatrix
    acc: ComplexMatrix
    b: ComplexMatrix

    def __call__(self, s: complex) -> ComplexMatrix:
        identity = np.eye(self.acc.shape[0], dtype=np.complex128)
        return self.d0 + self.c @ np.linalg.solve(s * identity - self.acc, self.b)

