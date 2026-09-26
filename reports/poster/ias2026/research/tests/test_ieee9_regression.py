"""GATE 1. The IEEE-9 GFL case, built independently from the stated equations."""

from __future__ import annotations

import pytest

from ibr_cycles.actions.interaction import split_interaction
from ibr_cycles.actions.portfolios import subsets
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee9_case import (
    ABSCISSA_TOLERANCE,
    ACTION_NAMES,
    EXPECTED_ABSCISSA,
    EXPECTED_CRITICAL_IMAG,
    PROBE_POINTS,
    SPECIFICATION_ABSCISSA,
    SPECIFICATION_TOLERANCE,
    build_case,
)
from ibr_cycles.reduction.schur import schur_determinant_residual


@pytest.fixture(scope="module")
def ieee9():
    return build_case()


def test_power_flow_is_converged_and_within_voltage_limits(ieee9):
    assert ieee9.power_flow.converged
    assert ieee9.power_flow.max_mismatch < 1e-10
    assert 0.95 < ieee9.power_flow.min_voltage
    assert ieee9.power_flow.max_voltage < 1.06


def test_every_portfolio_equilibrium_is_certified(ieee9):
    for members in subsets(ACTION_NAMES):
        equilibrium = ieee9.equilibrium(members)
        assert equilibrium.ok, (members, equilibrium.status)
        assert equilibrium.norm_f < 1e-9
        assert equilibrium.norm_g < 1e-9


@pytest.mark.parametrize("members", list(EXPECTED_ABSCISSA))
def test_frozen_spectral_abscissa(ieee9, members):
    spectrum = eigen_analysis(ieee9.state_matrix(members))
    assert spectrum.spectral_abscissa == pytest.approx(
        EXPECTED_ABSCISSA[members], abs=ABSCISSA_TOLERANCE
    )


@pytest.mark.parametrize("members", list(SPECIFICATION_ABSCISSA))
def test_independent_corroboration_of_the_specification_envelope(ieee9, members):
    """Not a reproduction: an agreement measurement against an outside model."""

    spectrum = eigen_analysis(ieee9.state_matrix(members))
    assert spectrum.spectral_abscissa == pytest.approx(
        SPECIFICATION_ABSCISSA[members], abs=SPECIFICATION_TOLERANCE
    )


def test_minimum_destabilizing_order_is_three(ieee9):
    for members in subsets(ACTION_NAMES):
        unstable = eigen_analysis(ieee9.state_matrix(members)).spectral_abscissa > 0.0
        assert unstable == (len(members) == 3)


def test_created_mode_frequency(ieee9):
    critical = eigen_analysis(ieee9.state_matrix(ACTION_NAMES)).critical
    assert abs(critical.imag) == pytest.approx(EXPECTED_CRITICAL_IMAG, abs=1e-5)


def test_operating_point_is_invariant_under_every_action(ieee9):
    """All three actions scale a quantity that vanishes at equilibrium."""

    assert ieee9.equilibrium_drift() < 1e-10


def test_reduced_matrix_is_additive_in_the_actions(ieee9):
    """The low-rank calculus needs this; it is measured, not assumed."""

    assert ieee9.additivity_residual() < 1e-12


def test_every_action_is_rank_one(ieee9):
    for factor in ieee9.factors():
        assert factor.numerical_rank == 1
        assert factor.reconstruction_error < 1e-12


def test_schur_identity_on_the_converter_partition(ieee9):
    residual = schur_determinant_residual(ieee9.a0, ieee9.partition, PROBE_POINTS)
    assert residual.max_relative_error < 1e-10


def test_determinant_lemma_and_factorization(ieee9):
    green = ieee9.green()
    assert green.lemma_residual(PROBE_POINTS) < 1e-10
    for s in PROBE_POINTS:
        assert split_interaction(green, s).factorization_error < 1e-10


def test_created_mode_is_a_zero_of_the_collective_factor(ieee9):
    """At the created eigenvalue the collective factor vanishes, not the individual."""

    green = ieee9.green()
    critical = eigen_analysis(ieee9.state_matrix(ACTION_NAMES)).critical.value
    split = split_interaction(green, critical)
    assert abs(split.full) < 1e-10
    assert abs(split.collective) < 1e-8
    assert abs(split.individual) > 1e-3
    assert split.min_singular_value < 1e-10
