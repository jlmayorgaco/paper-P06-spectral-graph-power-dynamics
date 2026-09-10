"""EXACT IDENTITY: matrix determinant lemma and the individual/collective split."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.actions.interaction import split_interaction
from ibr_cycles.models.toy5_case import ACTION_NAMES, PROBE_POINTS


def test_every_toy_action_is_exactly_rank_one(toy):
    for factor in toy.factors:
        assert factor.numerical_rank == 1
        assert factor.reconstruction_error < 1e-14
        assert factor.singular_values[1] < 1e-14 * factor.singular_values[0]


@pytest.mark.parametrize("members", [("A",), ("A", "B"), ACTION_NAMES])
def test_determinant_lemma(toy, members):
    green = toy.green(members)
    assert green.lemma_residual(PROBE_POINTS) < 1e-10


@pytest.mark.parametrize("members", [("A",), ("A", "B"), ACTION_NAMES])
def test_lemma_matches_direct_state_matrix(toy, members):
    green = toy.green(members)
    for s in PROBE_POINTS:
        direct = np.linalg.det(s * np.eye(5) - toy.state_matrix(members))
        base = np.linalg.det(s * np.eye(5) - toy.a0)
        assert abs(green.determinant_ratio(s) - direct / base) < 1e-10


def test_individual_times_collective_reproduces_full(toy):
    green = toy.green(ACTION_NAMES)
    for s in PROBE_POINTS:
        split = split_interaction(green, s)
        assert split.factorization_error < 1e-12


def test_individual_factor_is_the_product_of_single_action_ratios(toy):
    green = toy.green(ACTION_NAMES)
    for s in PROBE_POINTS:
        split = split_interaction(green, s)
        product = 1.0 + 0.0j
        for name in ACTION_NAMES:
            product *= toy.green((name,)).determinant_ratio(s)
        assert abs(split.individual - product) < 1e-12
