from spectral_ibr.application.ports.nep_solver import INepSolver
from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.entities.reduced_object import ReducedObject
from spectral_ibr.domain.value_objects.contour import Contour


def solve_nep_by_contour(
    solver: INepSolver,
    reduced: ReducedObject,
    contour: Contour,
) -> list[Pole]:
    return solver.solve(reduced, contour)

