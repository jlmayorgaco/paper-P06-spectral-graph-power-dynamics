from dataclasses import dataclass


@dataclass(frozen=True)
class PoleSensitivity:
    ds_dp: complex
    dzeta_dp: float


def compute_pole_sensitivity(*_: object, **__: object) -> PoleSensitivity:
    raise NotImplementedError("NEP pole sensitivity will be migrated here.")

