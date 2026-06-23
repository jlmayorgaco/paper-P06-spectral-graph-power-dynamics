from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PowerSystemCase:
    """Identity and provenance for a study case."""

    name: str
    source: Path | None = None
    base_mva: float | None = None
    tags: tuple[str, ...] = field(default_factory=tuple)

