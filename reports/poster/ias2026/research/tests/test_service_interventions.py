"""F8 service interventions are physical model changes with known limits.

Each switch must reduce to a documented configuration at its ends, and must not
move the operating point.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

CORE = (30, 33, 35, 37)


def _eig(case):
    return np.linalg.eigvals(case.system.A)


def _distance(a, b):
    cost = np.abs(a[:, None] - b[None, :])
    rows, cols = linear_sum_assignment(cost)
    return float(cost[rows, cols].max())


def test_avr_blend_zero_is_manual_excitation():
    blend = solve_case(ReplacementPlan.of({}), machine_services={"avr_blend": 0.0})
    manual = solve_case(ReplacementPlan.of({}), machine_services={"avr_manual": 1.0})
    assert np.allclose(blend.system.A, manual.system.A, atol=1e-7)


def test_flux_blend_zero_freezes_both_transient_emfs():
    case = solve_case(ReplacementPlan.of({}), machine_services={"flux_blend": 0.0})
    rows = [
        i for i, n in enumerate(case.system.labels) if n.startswith(("eq1_", "ed1_"))
    ]
    assert rows
    assert np.max(np.abs(case.system.A[rows])) == 0.0


def test_services_do_not_move_the_operating_point():
    plain = solve_case(ReplacementPlan.of({}))
    for services in (
        {"flux_blend": 0.3},
        {"avr_blend": 0.2},
        {"damping": 2.0},
        {"inertia": 0.5, "pss": 0.0},
    ):
        case = solve_case(ReplacementPlan.of({}), machine_services=services)
        assert np.allclose(case.equilibrium.z, plain.equilibrium.z, atol=1e-9)


def test_virtual_inertia_at_zero_gain_adds_one_decoupled_pole_per_converter():
    base = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    zero = solve_case(
        ReplacementPlan.of({b: 1.0 for b in CORE}),
        converter=ConverterParameters(
            inertia_emulation=True, h_virtual=0.0, inertia_filter=0.05
        ),
    )
    expected = np.concatenate([_eig(base), np.full(len(CORE), -1.0 / 0.05)])
    assert zero.n_states == base.n_states + len(CORE)
    assert _distance(_eig(zero), expected) < 1e-5


def test_virtual_inertia_changes_the_spectrum_when_active():
    zero = solve_case(
        ReplacementPlan.of({b: 1.0 for b in CORE}),
        converter=ConverterParameters(inertia_emulation=True),
    )
    active = solve_case(
        ReplacementPlan.of({b: 1.0 for b in CORE}),
        converter=ConverterParameters(inertia_emulation=True, h_virtual=4.0),
    )
    assert _distance(_eig(zero), _eig(active)) > 1e-3


def test_default_condenser_share_reproduces_the_frozen_condenser():
    kwargs = dict(condenser={b: 0.5 for b in CORE})
    old = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}, **kwargs))
    new = solve_case(
        ReplacementPlan.of(
            {b: 1.0 for b in CORE}, **kwargs, condenser_services={"q_share": 1.0}
        )
    )
    assert np.array_equal(old.system.A, new.system.A)


@pytest.mark.parametrize("share", [0.0, 0.5])
def test_condenser_reactive_share_keeps_the_network_solution(share):
    ref = solve_case(
        ReplacementPlan.of({b: 1.0 for b in CORE}, condenser={b: 0.5 for b in CORE})
    )
    case = solve_case(
        ReplacementPlan.of(
            {b: 1.0 for b in CORE},
            condenser={b: 0.5 for b in CORE},
            condenser_services={"q_share": share},
        )
    )
    # same network solution; differences are the equilibrium solver tolerance
    assert np.allclose(case.equilibrium.z, ref.equilibrium.z, atol=1e-6)
