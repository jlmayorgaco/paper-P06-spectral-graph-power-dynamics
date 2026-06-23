from spectral_ibr.domain.entities.linearized_model import LinearizedModel
from spectral_ibr.domain.entities.pole import Pole


class AndesEigensolver:
    def solve(self, model: LinearizedModel) -> list[Pole]:
        raise NotImplementedError("Wire ANDES eigenanalysis here.")

