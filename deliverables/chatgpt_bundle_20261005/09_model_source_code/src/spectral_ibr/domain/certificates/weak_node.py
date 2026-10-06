from dataclasses import dataclass


@dataclass(frozen=True)
class WeakNodeCertificate:
    node_id: str
    score: float
    margin: float

    @property
    def is_weak(self) -> bool:
        return self.score > self.margin

