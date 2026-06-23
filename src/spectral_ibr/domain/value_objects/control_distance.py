from dataclasses import dataclass


@dataclass(frozen=True)
class ControlDistance:
    """Control distance and passivity-defect observables."""

    dc: float
    phi_c: float
    passivity_defect: float

