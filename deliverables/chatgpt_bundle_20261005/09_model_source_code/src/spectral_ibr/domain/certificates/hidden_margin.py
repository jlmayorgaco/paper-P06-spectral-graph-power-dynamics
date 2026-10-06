from dataclasses import dataclass


@dataclass(frozen=True)
class HiddenMarginCertificate:
    mode_id: str
    observed_margin: float
    hidden_margin: float

    @property
    def is_hidden_limited(self) -> bool:
        return self.hidden_margin < self.observed_margin

