from dataclasses import dataclass, field
from enum import StrEnum


class ActionKind(StrEnum):
    LINE = "line"
    INERTIA = "inertia"
    DAMPING = "damping"
    CONTROL = "control"


@dataclass(frozen=True)
class Action:
    """Parametric intervention evaluated by a gate or certificate."""

    kind: ActionKind
    target: str
    magnitude: float
    metadata: dict[str, str | float] = field(default_factory=dict)

