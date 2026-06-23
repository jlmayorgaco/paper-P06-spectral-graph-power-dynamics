from typing import Protocol

from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.entities.reduced_object import ReducedObject
from spectral_ibr.domain.value_objects.contour import Contour


class INepSolver(Protocol):
    def solve(self, reduced: ReducedObject, contour: Contour) -> list[Pole]:
        """Solve a nonlinear eigenvalue problem inside a contour."""

