from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.entities.reduced_object import ReducedObject


class NewtonPolish:
    def polish(self, reduced: ReducedObject, pole: Pole) -> Pole:
        raise NotImplementedError("Implement non-circular NEP pole polishing here.")

