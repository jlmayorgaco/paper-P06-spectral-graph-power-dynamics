"""Frozen IEEE-9 GFL case: equilibrium, reduction, action factorization.

Two envelopes are stored side by side and they are not the same thing.

``EXPECTED_ABSCISSA`` is what *this* implementation produces; it is the
regression contract and is checked at ``ABSCISSA_TOLERANCE``.

``SPECIFICATION_ABSCISSA`` is the envelope quoted in the original research
specification. Its generating script is not in this repository and the textual
description does not pin the model to six decimals, so it cannot be reproduced
exactly. It is kept as an independent corroboration target at
``SPECIFICATION_TOLERANCE``: agreement at that level on all eight portfolios,
including the sign change at the triple, is evidence that the two independent
implementations describe the same system, and it is reported as such rather
than as a reproduction.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache, lru_cache

import numpy as np
from numpy.typing import NDArray

from ..actions.green import ActionGreen
from ..actions.low_rank import ActionFactor, factorize_state_update
from ..dynamics.equilibrium import Equilibrium, solve_equilibrium
from ..dynamics.linearize import central_difference_jacobians, reduce_index_one
from ..reduction.self_energy import StatePartition
from . import ieee9_gfl as model

ACTION_NAMES = ("A", "B", "C")

EXPECTED_ABSCISSA: dict[tuple[str, ...], float] = {
    (): -0.105890,
    ("A",): -0.027533,
    ("B",): -0.086124,
    ("C",): -0.112056,
    ("A", "B"): -0.008014,
    ("A", "C"): -0.034409,
    ("B", "C"): -0.073377,
    ("A", "B", "C"): +0.003217,
}
EXPECTED_BASE_IMAG = 8.578853
EXPECTED_CRITICAL_IMAG = 8.530945
ABSCISSA_TOLERANCE = 5e-6

SPECIFICATION_ABSCISSA: dict[tuple[str, ...], float] = {
    (): -0.105542,
    ("A",): -0.027190,
    ("B",): -0.085958,
    ("C",): -0.111529,
    ("A", "B"): -0.007848,
    ("A", "C"): -0.033884,
    ("B", "C"): -0.073325,
    ("A", "B", "C"): +0.003286,
}
SPECIFICATION_CRITICAL_IMAG = 8.531267
SPECIFICATION_TOLERANCE = 6e-4

#: Probe frequencies for the exact identities, away from every spectrum.
PROBE_POINTS: tuple[complex, ...] = (
    0.50 + 8.50j,
    -2.00 + 3.00j,
    1.00 - 12.00j,
    -15.00 + 20.00j,
)

#: Region boundary that isolates the portfolio modes from every subset mode.
RHP_BOUNDARY = -0.002
RHP_RADIUS = 30.0


@dataclass(frozen=True, eq=False)
class Ieee9Case:
    """The frozen operating point with everything derived from it.

    ``eq=False`` keeps identity hashing so the reduction cache can key on the
    case; the fields are numpy arrays and are not comparable by value.
    """

    parameters: model.Ieee9Parameters
    power_flow: model.PowerFlow
    x0: NDArray[np.float64]
    z0: NDArray[np.float64]

    @property
    def partition(self) -> StatePartition:
        retained = tuple(
            model.STATE_LABELS.index(name) for name in model.RETAINED_STATES
        )
        return StatePartition.from_retained(len(model.STATE_LABELS), retained)

    def equilibrium(self, members: tuple[str, ...]) -> Equilibrium:
        stressed = model.apply_actions(self.parameters, members)
        dae = model.Ieee9GflModel(parameters=stressed)
        return solve_equilibrium(dae, {}, self.x0, self.z0, tol=1e-9)

    def state_matrix(self, members: tuple[str, ...]) -> NDArray[np.float64]:
        return _reduced(self, members)

    @property
    def a0(self) -> NDArray[np.float64]:
        return self.state_matrix(())

    def factors(
        self, members: tuple[str, ...] = ACTION_NAMES
    ) -> tuple[ActionFactor, ...]:
        base = self.a0
        return tuple(
            factorize_state_update(
                name, self.state_matrix((name,)) - base, relative_tolerance=1e-6
            )
            for name in members
        )

    def green(self, members: tuple[str, ...] = ACTION_NAMES) -> ActionGreen:
        return ActionGreen(
            a0=self.a0.astype(np.complex128), factors=self.factors(members)
        )

    def additivity_residual(self) -> float:
        """How far the reduced matrix is from being additive in the actions.

        The three actions enter different rows of ``f`` and never enter ``g``,
        so the index-1 reduction stays affine in each of them and the portfolio
        update is the sum of the single-action updates. This measures it rather
        than assuming it, because the low-rank calculus is built on it.
        """

        base = self.a0
        predicted = base + sum(
            self.state_matrix((name,)) - base for name in ACTION_NAMES
        )
        actual = self.state_matrix(ACTION_NAMES)
        return float(
            np.abs(predicted - actual).max() / max(float(np.abs(actual).max()), 1e-300)
        )

    def equilibrium_drift(self) -> float:
        """Largest movement of the operating point under any single action.

        All three actions multiply a quantity that vanishes at the equilibrium
        (speed deviation, power error, ``v_q``), so the operating point does not
        move and the frozen-point and re-equilibrated readings of METHODS.md
        section 10 coincide for this action set. That is a property of these
        actions, not a general licence.
        """

        drift = 0.0
        for members in [(), ("A",), ("B",), ("C",), ACTION_NAMES]:
            equilibrium = self.equilibrium(members)
            drift = max(
                drift,
                float(np.abs(equilibrium.x - self.x0).max()),
                float(np.abs(equilibrium.z - self.z0).max()),
            )
        return drift


@cache
def _reduced(case: Ieee9Case, members: tuple[str, ...]) -> NDArray[np.float64]:
    stressed = model.apply_actions(case.parameters, members)
    dae = model.Ieee9GflModel(parameters=stressed)
    equilibrium = solve_equilibrium(dae, {}, case.x0, case.z0, tol=1e-9)
    if not equilibrium.ok:
        raise RuntimeError(f"equilibrium failed for {members}: {equilibrium.status}")
    jacobians = central_difference_jacobians(dae, equilibrium.x, equilibrium.z, {})
    return reduce_index_one(jacobians, model.STATE_LABELS).A


@lru_cache(maxsize=1)
def build_case() -> Ieee9Case:
    power_flow = model.solve_power_flow()
    if not power_flow.converged:
        raise RuntimeError("base power flow did not converge")
    parameters, x0, z0 = model.initialize(power_flow, model.Ieee9Parameters())
    return Ieee9Case(parameters=parameters, power_flow=power_flow, x0=x0, z0=z0)
