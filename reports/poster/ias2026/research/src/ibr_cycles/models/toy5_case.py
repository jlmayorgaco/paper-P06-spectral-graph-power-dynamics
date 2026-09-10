from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..actions.green import ActionGreen
from ..actions.low_rank import ActionFactor, factorize_state_update
from ..reduction.self_energy import StatePartition
from . import toy5

ACTION_NAMES = ("A", "B", "C")

#: Frozen regression envelope. Every number is Re(lambda_crit) of the portfolio.
EXPECTED_ABSCISSA: dict[tuple[str, ...], float] = {
    (): -0.199881,
    ("A",): -0.204175,
    ("B",): -0.225588,
    ("C",): -0.034298,
    ("A", "B"): -0.191010,
    ("A", "C"): -0.020165,
    ("B", "C"): -0.076228,
    ("A", "B", "C"): +0.046781,
}
EXPECTED_CRITICAL_IMAG = 2.212347
EXPECTED_K_CRITICAL = 1.351826
EXPECTED_K_REPAIR = 1.259557
ABSCISSA_TOLERANCE = 1e-6

#: Probe frequencies for the exact identities; chosen away from any spectrum.
PROBE_POINTS: tuple[complex, ...] = (
    0.30 + 1.10j,
    -0.20 + 1.94j,
    0.05 + 2.21j,
    2.00 - 0.70j,
    -1.30 + 0.00j,
)


@dataclass(frozen=True)
class ToyCase:
    """The frozen five-state regression case with its action factorization."""

    parameters: toy5.ToyParameters
    factors: tuple[ActionFactor, ...]

    @property
    def a0(self) -> np.ndarray:
        return toy5.state_matrix(self.parameters)

    @property
    def partition(self) -> StatePartition:
        return StatePartition.from_retained(5, toy5.RETAINED)

    def green(self, members: tuple[str, ...] = ACTION_NAMES) -> ActionGreen:
        chosen = tuple(f for f in self.factors if f.name in members)
        return ActionGreen(a0=self.a0.astype(np.complex128), factors=chosen)

    def state_matrix(self, members: tuple[str, ...]) -> np.ndarray:
        return toy5.state_matrix(toy5.apply_actions(self.parameters, members))


def build_case(parameters: toy5.ToyParameters | None = None) -> ToyCase:
    base = parameters or toy5.ToyParameters()
    factors = tuple(
        factorize_state_update(name, toy5.delta_a(base, name)) for name in ACTION_NAMES
    )
    return ToyCase(parameters=base, factors=factors)
