"""Mode provenance, and the contour dependence that makes it conditional."""

from __future__ import annotations

from ibr_cycles.actions.portfolios import subsets
from ibr_cycles.cycles.winding import (
    check_admissibility,
    circle,
    right_half_plane,
    spectra_union,
    winding_number,
)
from ibr_cycles.diagnosis.provenance import decompose_winding
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.toy5_case import ACTION_NAMES

CREATED_MODE = 0.046781 + 2.212347j


def _spectrum(toy, members):
    return eigen_analysis(toy.state_matrix(members)).values


def _admissibility(toy, contour):
    lower = spectra_union(
        [
            _spectrum(toy, m)
            for m in subsets(ACTION_NAMES, include_empty=False)
            if len(m) < 3
        ]
    )
    return check_admissibility(
        contour,
        portfolio_modes=_spectrum(toy, ACTION_NAMES),
        base_modes=_spectrum(toy, ()),
        lower_order_modes=lower,
    )


def test_right_half_plane_counts_every_unstable_eigenvalue(toy):
    """A conjugate pair contributes +2. The count is per eigenvalue."""

    contour = right_half_plane(boundary=-5e-3, radius=60.0, samples=8192)
    green = toy.green(ACTION_NAMES)
    full = winding_number(lambda s, g=green: g.determinant_ratio(s), contour)
    assert full.nearest_integer == 2
    assert full.converged


def test_no_proper_subset_creates_a_right_half_plane_mode(toy):
    contour = right_half_plane(boundary=-5e-3, radius=60.0, samples=8192)
    for members in subsets(ACTION_NAMES, include_empty=False):
        if len(members) == 3:
            continue
        green = toy.green(members)
        result = winding_number(lambda s, g=green: g.determinant_ratio(s), contour)
        assert result.nearest_integer == 0


def test_provenance_is_collective_on_an_admissible_region(toy):
    contour = circle(CREATED_MODE, 0.15, samples=2048)
    admissibility = _admissibility(toy, contour)
    assert admissibility.admissible
    provenance = decompose_winding(toy.green(ACTION_NAMES), contour, admissibility)
    assert provenance.consistent
    assert provenance.converged
    assert provenance.full.nearest_integer == 1
    assert provenance.individual.nearest_integer == 0
    assert provenance.collective.nearest_integer == 1
    assert provenance.verdict == "COLLECTIVE"


def test_provenance_is_refused_on_an_inadmissible_region(toy):
    """A larger circle swallows the single-action mode of A and flips the split.

    Both windings stay arithmetically correct, so the guard has to be the
    admissibility test, not the integer residual.
    """

    contour = circle(CREATED_MODE, 0.30, samples=2048)
    admissibility = _admissibility(toy, contour)
    assert not admissibility.admissible
    assert admissibility.enclosed_lower_order >= 1
    provenance = decompose_winding(toy.green(ACTION_NAMES), contour, admissibility)
    assert provenance.verdict == "INADMISSIBLE_REGION"
    assert provenance.full.nearest_integer == 1
    assert provenance.individual.nearest_integer == 1
    assert provenance.collective.nearest_integer == 0
