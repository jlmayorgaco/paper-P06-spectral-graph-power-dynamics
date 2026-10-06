from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.entities.reduced_object import ReducedObject
from spectral_ibr.domain.value_objects.contour import Contour


class BeynContourSolver:
    def solve(self, reduced: ReducedObject, contour: Contour) -> list[Pole]:
        raise NotImplementedError("Implement Beyn contour moments here.")

