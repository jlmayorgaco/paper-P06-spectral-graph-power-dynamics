from dataclasses import dataclass, field

from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.value_objects.mode_family import ModeFamily


@dataclass(frozen=True)
class Mode:
    """Classified pole with participation metadata."""

    pole: Pole
    family: ModeFamily
    pi_c: float
    participation: dict[str, float] = field(default_factory=dict)

