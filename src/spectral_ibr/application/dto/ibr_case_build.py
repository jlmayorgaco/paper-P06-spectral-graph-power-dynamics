from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class IbrReplacementGroup:
    """Declarative SG-to-IBR replacement group."""

    kind: str
    gens: tuple[int, ...]
    buses: tuple[int, ...]
    model: str
    with_pll: bool = False
    parameter_set: str | None = None
    source: str | None = None


@dataclass(frozen=True)
class IbrConversionProtocol:
    """Citable conversion protocol for building a processed IBR case."""

    name: str
    base_case: Path
    output_case: Path
    profile_name: str
    eligible_buses: tuple[int, ...]
    keep_synchronous: tuple[int, ...]
    replacements: tuple[IbrReplacementGroup, ...]
    preserve_power_flow: bool
    parameter_sets: dict[str, dict[str, dict[str, object]]] = field(default_factory=dict)
    citations: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class IbrBuildRequest:
    case_name: str
    profile_name: str
    protocol_path: Path
    output_dir: Path
    seed: int = 0


@dataclass(frozen=True)
class IbrBuildResult:
    case_name: str
    status: str
    output_dir: Path
    manifest_path: Path
    artifacts: dict[str, Path] = field(default_factory=dict)
    metrics: dict[str, float] = field(default_factory=dict)
    message: str | None = None
