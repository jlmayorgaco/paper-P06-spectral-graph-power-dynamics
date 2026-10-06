from dataclasses import dataclass
from math import pi

from spectral_ibr.domain.value_objects.damping_margin import damping_ratio


@dataclass(frozen=True)
class Pole:
    """Eigenvalue and derived modal quantities."""

    s: complex
    frequency_hz: float
    damping_ratio: float
    residual: float | None = None

    @classmethod
    def from_eigenvalue(cls, s: complex, residual: float | None = None) -> "Pole":
        return cls(
            s=s,
            frequency_hz=abs(s.imag) / (2.0 * pi),
            damping_ratio=damping_ratio(s),
            residual=residual,
        )

