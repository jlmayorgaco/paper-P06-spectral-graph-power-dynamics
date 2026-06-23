from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from spectral_ibr.domain.errors import InvalidLinearizedModelError

FloatMatrix = NDArray[np.float64]


@dataclass(frozen=True)
class LinearizedModel:
    """Partitioned linearized dynamics matrix."""

    aee: FloatMatrix
    aec: FloatMatrix
    ace: FloatMatrix
    acc: FloatMatrix
    electrical_states: tuple[str, ...] = ()
    control_states: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        ne = self.aee.shape[0]
        nc = self.acc.shape[0]
        if self.aee.shape != (ne, ne):
            raise InvalidLinearizedModelError("aee must be square.")
        if self.acc.shape != (nc, nc):
            raise InvalidLinearizedModelError("acc must be square.")
        if self.aec.shape != (ne, nc):
            raise InvalidLinearizedModelError("aec must have shape (ne, nc).")
        if self.ace.shape != (nc, ne):
            raise InvalidLinearizedModelError("ace must have shape (nc, ne).")

    @property
    def full_matrix(self) -> FloatMatrix:
        return np.block([[self.aee, self.aec], [self.ace, self.acc]])

