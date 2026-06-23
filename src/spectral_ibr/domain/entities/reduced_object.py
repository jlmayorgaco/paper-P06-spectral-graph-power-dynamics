from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class ReducedObject:
    """Reduced nonlinear eigenvalue object S_e(s) or T(s)."""

    name: str
    evaluate: Callable[[complex], ComplexMatrix]

    def __call__(self, s: complex) -> ComplexMatrix:
        return self.evaluate(s)

