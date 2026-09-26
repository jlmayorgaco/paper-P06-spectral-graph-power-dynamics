"""Low-rank factorizations are not unique; only the invariants may be reported."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.cycles.holonomy import gauge_invariance_residual, holonomy
from ibr_cycles.models.toy5_case import ACTION_NAMES

PROBE = 0.05 + 2.21j


def _random_bases(rng, ranks):
    bases = []
    for rank in ranks:
        while True:
            candidate = rng.normal(size=(rank, rank)) + 1j * rng.normal(
                size=(rank, rank)
            )
            if np.linalg.cond(candidate) < 1e3:
                bases.append(candidate)
                break
    return bases


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_cycle_invariants_survive_random_basis_changes(toy, seed):
    rng = np.random.default_rng(seed)
    green = toy.green(ACTION_NAMES)
    bases = _random_bases(rng, [f.rank for f in green.factors])
    for cycle in [("A", "B"), ("A", "C"), ("B", "C"), ("A", "B", "C")]:
        drift = gauge_invariance_residual(green, cycle, PROBE, bases)
        assert drift["trace"] < 1e-9
        assert drift["determinant"] < 1e-9
        assert drift["eigenvalues"] < 1e-9


def test_holonomy_is_invariant_under_cyclic_rotation(toy):
    """THEOREM: the trace of a cyclic product does not depend on the base point."""

    green = toy.green(ACTION_NAMES)
    reference = holonomy(green, ("A", "B", "C"), PROBE).trace
    for rotation in [("B", "C", "A"), ("C", "A", "B")]:
        assert holonomy(green, rotation, PROBE).trace == pytest.approx(
            reference, rel=1e-12
        )


def test_toy_orientations_coincide(toy):
    """NUMERICAL OBSERVATION, this case only: H_ABC and H_ACB have equal traces.

    There is no general theorem forcing this, and the action Green operator here
    is demonstrably non-symmetric, so the coincidence is a property of the toy
    geometry. It is recorded rather than generalized, and the IEEE-9 and IEEE-39
    cases must be checked separately before any orientation claim is made.
    """

    green = toy.green(ACTION_NAMES)
    k = green.k(PROBE)
    assert np.abs(k - k.T).max() > 1e-3
    forward = holonomy(green, ("A", "B", "C"), PROBE).trace
    backward = holonomy(green, ("A", "C", "B"), PROBE).trace
    assert forward == pytest.approx(backward, rel=1e-10)
