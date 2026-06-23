from abc import ABC, abstractmethod
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True)
class GateRunConfig:
    """Runtime inputs shared by reproducibility gates."""

    name: str
    case_name: str | None = None
    profile_name: str | None = None
    config_path: Path | None = None
    output_root: Path = Path("outputs")
    seed: int = 0
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class GateResult:
    name: str
    verdict: str
    artifacts: dict[str, str] = field(default_factory=dict)
    metrics: dict[str, float] = field(default_factory=dict)
    message: str | None = None


class Gate(Protocol):
    name: str

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        """Execute the gate and return a reproducibility verdict."""


class BaseGate(ABC):
    """Template for gates that emit versioned artifacts."""

    name: str

    def setup(self, config: GateRunConfig) -> None:
        """Prepare directories or validate config before the gate runs."""

    @abstractmethod
    def run(self, config: GateRunConfig | None = None) -> GateResult:
        """Execute the gate."""

    def verdict(self, result: GateResult) -> str:
        return result.verdict

    def artifacts(self, result: GateResult) -> dict[str, str]:
        return result.artifacts
