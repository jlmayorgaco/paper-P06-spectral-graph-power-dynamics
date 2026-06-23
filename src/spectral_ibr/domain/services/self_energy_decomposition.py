from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class SelfEnergyTerms:
    direct: NDArray[np.complex128]
    coupling: NDArray[np.complex128]
    control: NDArray[np.complex128]


def decompose_self_energy(*_: object, **__: object) -> SelfEnergyTerms:
    raise NotImplementedError("Self-energy numerator decomposition belongs here.")

