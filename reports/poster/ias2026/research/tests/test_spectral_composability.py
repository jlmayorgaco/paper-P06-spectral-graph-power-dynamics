"""Spectral composability order: invariance, and its independence from interaction order.

The counterexamples here are the formal content of
``theory/spectral_composability_order.md`` section 4. They are kept as tests so
that the two notions can never quietly be conflated again in code.
"""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.diagnosis.composability import (
    NO_FINITE_ORDER,
    BandRightHalfPlane,
    composability_order,
    interaction_order,
    mobius_terms,
    spectral_count,
    witnesses,
)

BAND = BandRightHalfPlane(low_hz=0.3, high_hz=1.5)


def _pair(real, frequency_hz):
    omega = 2.0 * np.pi * frequency_hz
    return np.array([complex(real, omega), complex(real, -omega)])


# --------------------------------------------------------------- the region --

def test_region_counts_one_representative_per_conjugate_pair():
    assert spectral_count(_pair(+0.1, 0.6), BAND) == 1
    assert spectral_count(_pair(-0.1, 0.6), BAND) == 0


def test_region_excludes_modes_outside_the_band():
    assert spectral_count(_pair(+0.1, 0.2), BAND) == 0
    assert spectral_count(_pair(+0.1, 2.0), BAND) == 0


def test_region_excludes_the_reference_zero_modes():
    assert spectral_count(np.array([0.0 + 0.0j, 1e-9 + 0.0j]), BAND) == 0


def test_distance_to_boundary_changes_sign_exactly_at_a_crossing():
    assert BAND.distance_to_boundary(_pair(-0.05, 0.6)) < 0.0
    assert BAND.distance_to_boundary(_pair(+0.05, 0.6)) > 0.0


# ------------------------------------------------------------- invariance ----

def test_count_is_invariant_under_similarity():
    rng = np.random.default_rng(0)
    block = np.array([[-0.1, 4.0], [-4.0, -0.1]])
    matrix = np.zeros((4, 4))
    matrix[:2, :2] = block
    matrix[2:, 2:] = np.array([[+0.05, 4.4], [-4.4, +0.05]])
    transform = rng.normal(size=(4, 4))
    while abs(np.linalg.det(transform)) < 1e-3:
        transform = rng.normal(size=(4, 4))
    similar = np.linalg.solve(transform, matrix @ transform)
    assert spectral_count(np.linalg.eigvals(matrix), BAND) == spectral_count(
        np.linalg.eigvals(similar), BAND
    )


def test_count_is_invariant_under_eigenvalue_relabelling():
    values = np.concatenate([_pair(+0.1, 0.6), _pair(-0.2, 0.9)])
    rng = np.random.default_rng(3)
    assert spectral_count(values, BAND) == spectral_count(rng.permutation(values), BAND)


def test_order_needs_the_baseline_present():
    with pytest.raises(ValueError):
        composability_order({frozenset({1}): 1})


# ---------------------------------------------------- the order itself -------

def test_order_is_the_smallest_set_that_changes_the_count():
    counts = {
        frozenset(): 0,
        frozenset({1}): 0, frozenset({2}): 0, frozenset({3}): 0,
        frozenset({1, 2}): 0, frozenset({1, 3}): 0, frozenset({2, 3}): 0,
        frozenset({1, 2, 3}): 1,
    }
    assert composability_order(counts) == 3
    assert witnesses(counts) == (frozenset({1, 2, 3}),)


def test_order_is_infinite_when_nothing_changes():
    counts = {frozenset(): 2, frozenset({1}): 2, frozenset({1, 2}): 2}
    assert composability_order(counts) == NO_FINITE_ORDER
    assert witnesses(counts) == ()


def test_a_decrease_in_the_count_also_counts_as_a_change():
    """kappa asks whether the count DIFFERS, not whether it grows."""

    counts = {frozenset(): 1, frozenset({1}): 1, frozenset({1, 2}): 0}
    assert composability_order(counts) == 2


# ------------- composability order is NOT the interaction order --------------

def _threshold_counts(values, level=0.0):
    return {s: int(v > level) for s, v in values.items()}


def test_kappa_three_with_no_interaction_at_any_order():
    """A purely ADDITIVE set function can still have kappa = 3.

    Every Moebius term of order two and three vanishes, so there is no joint
    effect whatsoever, yet no single and no pair crosses while the triple does.
    kappa is a THRESHOLD notion, not an interaction notion.
    """

    values = {
        frozenset(): -3.0,
        frozenset({1}): -1.8, frozenset({2}): -1.8, frozenset({3}): -1.8,
        frozenset({1, 2}): -0.6, frozenset({1, 3}): -0.6, frozenset({2, 3}): -0.6,
        frozenset({1, 2, 3}): +0.6,
    }
    terms = mobius_terms(values)
    for subset, term in terms.items():
        if len(subset) == 1:
            assert term == pytest.approx(1.2)
        elif len(subset) >= 2:
            assert term == pytest.approx(0.0, abs=1e-12)
    assert interaction_order(values, tolerance=1e-9) == 1
    assert composability_order(_threshold_counts(values)) == 3


def test_kappa_one_with_a_large_third_order_interaction():
    """And the converse: a big joint term with kappa = 1."""

    values = {
        frozenset(): -1.0,
        frozenset({1}): +1.0, frozenset({2}): -1.0, frozenset({3}): -1.0,
        frozenset({1, 2}): +1.0, frozenset({1, 3}): +1.0, frozenset({2, 3}): -1.0,
        frozenset({1, 2, 3}): +5.0,
    }
    assert mobius_terms(values)[frozenset({1, 2, 3})] == pytest.approx(4.0)
    assert interaction_order(values, tolerance=1e-9) == 3
    assert composability_order(_threshold_counts(values)) == 1


def test_the_two_notions_agree_on_the_flagship_pattern():
    """They coincide on the pattern Track A actually exhibits, which is why the
    conflation went unnoticed: all proper subsets safe AND a large top term."""

    values = {
        frozenset(): -0.126,
        frozenset({1}): -0.139, frozenset({2}): -0.159,
        frozenset({3}): -0.159, frozenset({4}): -0.130,
        frozenset({1, 2}): -0.206, frozenset({1, 3}): -0.194,
        frozenset({1, 4}): -0.153, frozenset({2, 3}): -0.233,
        frozenset({2, 4}): -0.172, frozenset({3, 4}): -0.171,
        frozenset({1, 2, 3}): -0.163, frozenset({1, 2, 4}): -0.277,
        frozenset({1, 3, 4}): -0.238, frozenset({2, 3, 4}): -0.280,
        frozenset({1, 2, 3, 4}): +0.145,
    }
    assert composability_order(_threshold_counts(values)) == 4
    assert interaction_order(values, tolerance=1e-6) == 4
    assert mobius_terms(values)[frozenset({1, 2, 3, 4})] == pytest.approx(0.433, abs=5e-3)
