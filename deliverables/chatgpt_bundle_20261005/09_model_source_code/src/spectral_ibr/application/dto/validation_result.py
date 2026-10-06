from dataclasses import dataclass


@dataclass(frozen=True)
class ValidationResult:
    case_name: str
    max_frequency_error_hz: float
    max_damping_error: float
    family_agreement: float
    passed: bool

