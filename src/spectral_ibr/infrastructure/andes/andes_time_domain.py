from spectral_ibr.domain.entities.power_system_case import PowerSystemCase


class AndesTimeDomainSimulator:
    def simulate(self, case: PowerSystemCase, horizon_s: float) -> dict[str, object]:
        raise NotImplementedError(f"Wire ANDES simulation for {case.name}, {horizon_s}s.")

