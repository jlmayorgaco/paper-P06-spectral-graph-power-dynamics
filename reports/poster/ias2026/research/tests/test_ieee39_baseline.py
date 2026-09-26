"""PHASE A eligibility and the exactness of the matched-policy frozen point."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

BASE_ZETA_MIN = 0.026124
BASE_ABSCISSA = -0.126478


@pytest.fixture(scope="module")
def base():
    return solve_case(ReplacementPlan.of({}))


def test_base_dimensions(base):
    assert base.n_states == 70
    assert base.dae.n_z == 78


def test_base_equilibrium_and_conditioning(base):
    assert base.equilibrium.ok
    assert base.equilibrium.norm_f < 1e-9
    assert base.equilibrium.norm_g < 1e-9
    assert base.gz_condition < 1e6


def test_base_is_eligible(base):
    verdict = assess(eigen_analysis(base.system.A), base.system.labels)
    assert not verdict.unstable
    assert verdict.zeta_min == pytest.approx(BASE_ZETA_MIN, abs=1e-5)
    assert verdict.spectral_abscissa == pytest.approx(BASE_ABSCISSA, abs=1e-5)


def test_reference_modes_are_the_free_angle_and_frequency(base):
    verdict = assess(eigen_analysis(base.system.A), base.system.labels)
    assert verdict.reference_modes == 2
    assert verdict.reference_participation > 0.99
    assert verdict.reference_max_abs < 1e-3


def test_matched_policy_leaves_the_operating_point_exactly_fixed(base):
    """The replacement keeps P and Q at the bus, so the AC solution cannot move.

    Two statements, of different strength. The power-flow solution is identical
    bit for bit, because the scheduled injections are unchanged: this is exact,
    not an approximation, and it is what makes the frozen-operator campaign
    legitimate. The DAE equilibrium then agrees only to the residual tolerance of
    its own solver, which is a property of the solver and not of the physics.
    """

    for members in [(30,), (33, 35), (30, 33, 35, 38)]:
        case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
        moved = np.abs(
            case.dae.power_flow.voltages - base.dae.power_flow.voltages
        ).max()
        assert moved == 0.0, members
        assert np.abs(case.equilibrium.z - base.equilibrium.z).max() < 1e-6


def test_every_single_replacement_is_feasible_and_stable():
    for bus in (30, 31, 32, 33, 34, 35, 36, 37, 38):
        case = solve_case(ReplacementPlan.of({bus: 1.0}))
        verdict = assess(eigen_analysis(case.system.A), case.system.labels)
        assert not verdict.unstable, bus
