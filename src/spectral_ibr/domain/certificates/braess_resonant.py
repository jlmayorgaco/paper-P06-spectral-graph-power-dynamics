from dataclasses import dataclass


@dataclass(frozen=True)
class BraessResonantCertificate:
    action_id: str
    sensitivity: float
    resonance_gap: float
    threshold: float

    @property
    def predicts_adverse_shift(self) -> bool:
        return self.sensitivity > 0.0 and self.resonance_gap <= self.threshold

