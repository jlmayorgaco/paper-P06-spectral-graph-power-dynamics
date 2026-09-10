"""Finite-amplitude interaction orders, and whether order three actually matters."""

from __future__ import annotations

import pytest

from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.toy5_case import ACTION_NAMES, EXPECTED_ABSCISSA


@pytest.fixture(scope="module")
def abscissa_decomposition(toy):
    def functional(members):
        return eigen_analysis(toy.state_matrix(tuple(members))).spectral_abscissa

    return decompose(functional, ACTION_NAMES)


def test_mobius_reconstruction_is_exact(abscissa_decomposition):
    assert abscissa_decomposition.reconstruction_error(ACTION_NAMES) < 1e-12


def test_pairwise_truncation_predicts_a_stable_portfolio(abscissa_decomposition):
    """The triple term is not an accumulation of lower orders.

    NUMERICAL OBSERVATION, frozen amplitudes only: truncating at order two puts
    the portfolio at -0.0232, comfortably stable. The irreducible third-order
    term of +0.0700 is what crosses the imaginary axis.
    """

    order_one = abscissa_decomposition.truncated(ACTION_NAMES, 1)
    order_two = abscissa_decomposition.truncated(ACTION_NAMES, 2)
    order_three = abscissa_decomposition.truncated(ACTION_NAMES, 3)
    assert order_one < 0.0
    assert order_two < 0.0
    assert order_three > 0.0
    assert order_three == pytest.approx(EXPECTED_ABSCISSA[ACTION_NAMES], abs=1e-6)


def test_third_order_term_magnitude(abscissa_decomposition):
    triple = abscissa_decomposition.terms[frozenset(ACTION_NAMES)]
    assert triple == pytest.approx(0.070005, abs=1e-5)
