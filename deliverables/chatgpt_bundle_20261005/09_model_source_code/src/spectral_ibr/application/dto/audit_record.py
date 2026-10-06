from dataclasses import dataclass, field


@dataclass(frozen=True)
class AuditRecord:
    gate: str
    case_name: str
    verdict: str
    metrics: dict[str, float] = field(default_factory=dict)

