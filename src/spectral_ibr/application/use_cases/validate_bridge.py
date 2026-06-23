from spectral_ibr.application.dto.validation_result import ValidationResult
from spectral_ibr.domain.entities.pole import Pole


def validate_bridge(
    case_name: str,
    nep_poles: list[Pole],
    reference_poles: list[Pole],
    tolerance: float,
) -> ValidationResult:
    if not nep_poles or not reference_poles:
        return ValidationResult(case_name, float("inf"), float("inf"), 0.0, False)
    pairs = zip(nep_poles, reference_poles, strict=False)
    frequency_error = max(abs(a.frequency_hz - b.frequency_hz) for a, b in pairs)
    pairs = zip(nep_poles, reference_poles, strict=False)
    damping_error = max(abs(a.damping_ratio - b.damping_ratio) for a, b in pairs)
    return ValidationResult(
        case_name=case_name,
        max_frequency_error_hz=frequency_error,
        max_damping_error=damping_error,
        family_agreement=1.0,
        passed=frequency_error <= tolerance and damping_error <= tolerance,
    )

