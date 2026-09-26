"""Safeguard A: the F7 reactive-policy coordinate must be regular at gain zero.

The naive coordinate scales the PI voltage-regulator gains to zero. That leaves
one decoupled integrator per converter, i.e. an artificial zero eigenvalue and a
continuum of equilibria, so gain zero cannot sit inside Theta_reg. The leaky
regulator used by F7 has no such point: at gain zero it is exactly the matched
fixed-Q system plus one decoupled stable pole per converter.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

CORE = (30, 33, 35, 37)
LEAK = 0.05


def _spectrum(converter=None, members=CORE):
    case = solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}), converter=converter
    )
    return case, np.linalg.eigvals(case.system.A)


def _matched_distance(a, b):
    cost = np.abs(a[:, None] - b[None, :])
    rows, cols = linear_sum_assignment(cost)
    return float(cost[rows, cols].max())


@pytest.fixture(scope="module")
def matched():
    return _spectrum()


def test_defaults_reproduce_the_plain_pi_regulator_bit_for_bit():
    left, _ = _spectrum(ConverterParameters(voltage_control=True))
    right, _ = _spectrum(
        ConverterParameters(voltage_control=True, voltage_gain=1.0, voltage_leak=0.0)
    )
    assert np.array_equal(left.system.A, right.system.A)


def test_naive_gain_zero_creates_one_zero_eigenvalue_per_converter(matched):
    _, values = _spectrum(ConverterParameters(voltage_control=True, voltage_gain=0.0))
    _, reference = matched
    extra_zeros = int(np.sum(np.abs(values) < 1e-9)) - int(
        np.sum(np.abs(reference) < 1e-9)
    )
    assert extra_zeros == len(CORE)


def test_leaky_gain_zero_is_matched_fixed_q_plus_decoupled_stable_poles(matched):
    case, values = _spectrum(
        ConverterParameters(voltage_control=True, voltage_gain=0.0, voltage_leak=LEAK)
    )
    matched_case, reference = matched
    assert case.n_states == matched_case.n_states + len(CORE)
    expected = np.concatenate([reference, np.full(len(CORE), -LEAK)])
    # central-difference Jacobians limit agreement to about 1e-6
    assert _matched_distance(values, expected) < 1e-5
    assert np.sum(np.abs(values) < 1e-9) == 0


def test_equilibrium_and_algebraic_block_do_not_depend_on_the_policy_gain():
    cases = [
        solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE}),
            converter=ConverterParameters(
                voltage_control=True, voltage_gain=g, voltage_leak=LEAK
            ),
        )
        for g in (0.0, 0.5, 1.5)
    ]
    for case in cases[1:]:
        assert np.allclose(case.equilibrium.x, cases[0].equilibrium.x, atol=1e-10)
        assert np.allclose(case.equilibrium.z, cases[0].equilibrium.z, atol=1e-10)
        assert case.gz_condition == pytest.approx(cases[0].gz_condition, rel=1e-8)


def test_reduced_matrix_is_affine_in_the_policy_gain():
    def matrix(g):
        return solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE}),
            converter=ConverterParameters(
                voltage_control=True, voltage_gain=g, voltage_leak=LEAK
            ),
        ).system.A

    a0, a1, ah = matrix(0.0), matrix(1.0), matrix(0.37)
    interpolated = (1.0 - 0.37) * a0 + 0.37 * a1
    assert np.linalg.norm(ah - interpolated) / np.linalg.norm(ah) < 1e-7


def test_bus_scaling_moves_one_machine_only():
    base = solve_case(ReplacementPlan.of({}))
    moved = solve_case(ReplacementPlan.of({}), machine_bus_scaling={34: {"ta": 2.0}})
    # rows other than the scaled one differ only by central-difference noise
    changed = np.flatnonzero(
        np.any(np.abs(base.system.A - moved.system.A) > 1e-6, axis=1)
    )
    labels = {base.system.labels[i] for i in changed}
    assert labels == {"efd_sg34"}
