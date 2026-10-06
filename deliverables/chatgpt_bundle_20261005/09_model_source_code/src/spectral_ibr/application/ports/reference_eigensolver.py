from typing import Protocol

from spectral_ibr.domain.entities.linearized_model import LinearizedModel
from spectral_ibr.domain.entities.pole import Pole


class IReferenceEigensolver(Protocol):
    def solve(self, model: LinearizedModel) -> list[Pole]:
        """Return reference eigenvalues from a full-order model."""

