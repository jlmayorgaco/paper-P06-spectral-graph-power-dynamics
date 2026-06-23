from dataclasses import dataclass


@dataclass(frozen=True)
class WeakLinkCertificate:
    link_id: str
    score: float
    margin: float

    @property
    def is_weak(self) -> bool:
        return self.score > self.margin

