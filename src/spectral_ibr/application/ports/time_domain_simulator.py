from typing import Protocol

from spectral_ibr.domain.entities.power_system_case import PowerSystemCase


class ITimeDomainSimulator(Protocol):
    def simulate(self, case: PowerSystemCase, horizon_s: float) -> dict[str, object]:
        """Run a time-domain simulation."""

