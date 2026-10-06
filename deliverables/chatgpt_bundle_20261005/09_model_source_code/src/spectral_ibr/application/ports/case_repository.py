from typing import Protocol

from spectral_ibr.domain.entities.linearized_model import LinearizedModel
from spectral_ibr.domain.entities.power_system_case import PowerSystemCase


class ICaseRepository(Protocol):
    def load(self, case: PowerSystemCase) -> LinearizedModel:
        """Load and partition a case."""

