"""Connected spectral cumulants must match the closed-form cycle identities."""

from __future__ import annotations

import pytest

from ibr_cycles.cycles.cumulants import (
    convergence_study,
    second_cumulant,
    third_cumulant,
)
from ibr_cycles.models.toy5_case import ACTION_NAMES

PROBE = 0.05 + 2.21j


@pytest.mark.parametrize("axes", [(0, 1), (0, 2), (1, 2)])
def test_second_cumulant_matches_finite_differences(toy, axes):
    green = toy.green(ACTION_NAMES)
    study = convergence_study(green, axes, PROBE)
    assert study.best_error < 1e-6


def test_third_cumulant_matches_finite_differences(toy):
    green = toy.green(ACTION_NAMES)
    study = convergence_study(green, (0, 1, 2), PROBE)
    assert study.best_error < 1e-5


def test_second_cumulant_is_symmetric(toy):
    green = toy.green(ACTION_NAMES)
    assert second_cumulant(green, 0, 1, PROBE) == pytest.approx(
        second_cumulant(green, 1, 0, PROBE), rel=1e-12
    )


def test_third_cumulant_is_permutation_symmetric(toy):
    green = toy.green(ACTION_NAMES)
    reference = third_cumulant(green, 0, 1, 2, PROBE)
    for a, b, c in [(1, 2, 0), (2, 0, 1), (0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        assert third_cumulant(green, a, b, c, PROBE) == pytest.approx(
            reference, rel=1e-10
        )


def test_second_cumulant_equals_minus_the_pair_holonomy(toy):
    """EXACT IDENTITY linking the cumulant expansion to closed action cycles."""

    from ibr_cycles.cycles.holonomy import holonomy

    green = toy.green(ACTION_NAMES)
    for a, b, names in [(0, 1, ("A", "B")), (0, 2, ("A", "C")), (1, 2, ("B", "C"))]:
        loop = holonomy(green, names, PROBE)
        assert second_cumulant(green, a, b, PROBE) == pytest.approx(
            -loop.trace, rel=1e-12
        )


def test_third_cumulant_equals_the_sum_of_both_triple_holonomies(toy):
    """EXACT IDENTITY: the third cumulant is trace(H_abc) + trace(H_acb)."""

    from ibr_cycles.cycles.holonomy import holonomy

    green = toy.green(ACTION_NAMES)
    forward = holonomy(green, ("A", "B", "C"), PROBE).trace
    backward = holonomy(green, ("A", "C", "B"), PROBE).trace
    assert third_cumulant(green, 0, 1, 2, PROBE) == pytest.approx(
        forward + backward, rel=1e-12
    )
