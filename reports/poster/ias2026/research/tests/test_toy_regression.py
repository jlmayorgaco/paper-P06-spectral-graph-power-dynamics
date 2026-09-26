"""GATE 0. The five-state case is the contract every later model inherits."""

from __future__ import annotations

from dataclasses import replace

import pytest

from ibr_cycles.actions.portfolios import subsets
from ibr_cycles.diagnosis.repair import scalar_repair
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models import toy5
from ibr_cycles.models.toy5_case import (
    ABSCISSA_TOLERANCE,
    ACTION_NAMES,
    EXPECTED_ABSCISSA,
    EXPECTED_CRITICAL_IMAG,
    EXPECTED_K_CRITICAL,
    EXPECTED_K_REPAIR,
)


@pytest.mark.parametrize("members", list(EXPECTED_ABSCISSA))
def test_portfolio_spectral_abscissa(toy, members):
    spectrum = eigen_analysis(toy.state_matrix(members))
    assert spectrum.spectral_abscissa == pytest.approx(
        EXPECTED_ABSCISSA[members], abs=ABSCISSA_TOLERANCE
    )


def test_full_portfolio_frequency(toy):
    spectrum = eigen_analysis(toy.state_matrix(ACTION_NAMES))
    assert abs(spectrum.critical.imag) == pytest.approx(
        EXPECTED_CRITICAL_IMAG, abs=ABSCISSA_TOLERANCE
    )


def test_minimum_destabilizing_order_is_three(toy):
    for members in subsets(ACTION_NAMES):
        spectrum = eigen_analysis(toy.state_matrix(members))
        unstable = spectrum.spectral_abscissa > 0.0
        assert unstable == (len(members) == 3)


def test_controller_surgery_endpoints(toy):
    def abscissa(k: float) -> float:
        stressed = toy5.apply_actions(toy.parameters, ("A", "B"))
        return float(
            eigen_analysis(toy5.state_matrix(replace(stressed, k=k))).spectral_abscissa
        )

    critical, repaired = scalar_repair(abscissa, (0.8, 1.57), target=-0.02)
    assert critical == pytest.approx(EXPECTED_K_CRITICAL, abs=1e-6)
    assert repaired == pytest.approx(EXPECTED_K_REPAIR, abs=1e-6)
