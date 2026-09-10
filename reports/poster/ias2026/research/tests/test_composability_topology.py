"""F2D: composability regions need not be monotone, nested or connected.

Analytic one-mode LTI families with additive actions,
    A(S; theta) = sigma(S; theta) I + omega0 J,  sigma = -1 + sum_{a in S} c_a(theta).
See theory/F2D_topology_of_composability_regions.md.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy import ndimage

from ibr_cycles.diagnosis.composability import (
    BandRightHalfPlane,
    composability_order,
    hypergraph_label,
    incompatibility_hypergraph,
    spectral_count,
    upward_closure_defect,
)

BAND = BandRightHalfPlane(low_hz=0.3, high_hz=1.5)
OMEGA0 = 2.0 * np.pi * 0.6
J = np.array([[0.0, 1.0], [-1.0, 0.0]])


def power_set(ground):
    return [
        frozenset(c) for k in range(len(ground) + 1) for c in combinations(ground, k)
    ]


def counts_for(c: dict, ground):
    out = {}
    for s in power_set(ground):
        sigma = -1.0 + sum(c[a] for a in s)
        out[s] = spectral_count(np.linalg.eigvals(sigma * np.eye(2) + OMEGA0 * J), BAND)
    return out


def test_kappa_falls_then_rises_along_a_transversal_path():
    lo, hi = (1 - 1 / np.sqrt(3)) / 2, (1 + 1 / np.sqrt(3)) / 2
    seen = []
    for theta in np.linspace(0.0, 1.0, 401):
        if min(abs(theta - lo), abs(theta - hi)) < 1e-6:
            continue
        counts = counts_for({1: 0.6 + 2.4 * theta * (1 - theta), 2: 0.6}, (1, 2))
        kappa = composability_order(counts)
        expected = 1 if lo < theta < hi else 2
        assert kappa == expected
        seen.append(kappa)
    runs = [k for i, k in enumerate(seen) if i == 0 or k != seen[i - 1]]
    assert runs == [2, 1, 2]
    # both crossings are transversal simple crossings
    assert abs(2.4 * (1 - 2 * lo)) > 1.0 and abs(2.4 * (1 - 2 * hi)) > 1.0


def test_single_action_infinity_one_infinity():
    seen = []
    for theta in np.linspace(0.0, 1.0, 201):
        counts = counts_for({1: 6 * theta * (1 - theta)}, (1,))
        seen.append(composability_order(counts))
    runs = [k for i, k in enumerate(seen) if i == 0 or k != seen[i - 1]]
    assert runs == [-1, 1, -1]


def test_order_one_region_has_two_islands_and_order_two_region_two_annuli():
    x = np.linspace(-2.0, 2.0, 321)
    xx, yy = np.meshgrid(x, x, indexing="ij")
    c1 = 1.3 - 4.0 * ((xx**2 - 1.0) ** 2 + yy**2)
    kappa = np.empty_like(c1, dtype=int)
    for idx in np.ndindex(c1.shape):
        kappa[idx] = composability_order(counts_for({1: c1[idx], 2: 0.6}, (1, 2)))
    islands, n_islands = ndimage.label(kappa == 1)
    assert n_islands == 2
    order_two, n_two = ndimage.label(kappa == 2)
    assert n_two == 2
    # holes: components of the complement of P_2 that do not touch the border
    complement, n_c = ndimage.label(kappa != 2)
    border = set(
        np.unique(
            np.concatenate(
                [complement[0], complement[-1], complement[:, 0], complement[:, -1]]
            )
        )
    ) - {0}
    assert n_c - len(border) == 2


def test_hypergraph_moves_while_kappa_stands_still():
    labels = []
    for theta in np.linspace(0.0, 1.0, 81):
        if min(abs(theta - 0.25), abs(theta - 0.5)) < 1e-9:
            continue
        counts = counts_for(
            {1: 0.6, 2: 0.6 - 0.4 * theta, 3: 0.3 + 0.4 * theta}, (1, 2, 3)
        )
        assert composability_order(counts) == 2
        labels.append(hypergraph_label(incompatibility_hypergraph(counts)))
    runs = [lab for i, lab in enumerate(labels) if i == 0 or lab != labels[i - 1]]
    assert runs == ["1+2", "1+2|1+3", "1+3"]


def test_one_stabilizing_action_breaks_heredity_without_interaction():
    counts = counts_for({1: 0.6, 2: 0.6, 3: -0.5}, (1, 2, 3))
    assert counts[frozenset({1, 2})] == 1
    assert counts[frozenset({1, 2, 3})] == 0
    assert frozenset({1, 2, 3}) in upward_closure_defect(counts)


def test_monotone_destabilisation_forces_monotone_kappa():
    rng = np.random.default_rng(7)
    for _ in range(50):
        slopes = rng.uniform(0.0, 1.0, 3)
        offsets = rng.uniform(-0.2, 0.6, 3)
        seen = []
        for theta in np.linspace(0.0, 1.0, 101):
            c = {a + 1: offsets[a] + slopes[a] * theta for a in range(3)}
            k = composability_order(counts_for(c, (1, 2, 3)))
            seen.append(np.inf if k == -1 else k)
        assert all(b <= a for a, b in zip(seen, seen[1:], strict=False))


def test_fold_normal_form_of_a_tongue_tip():
    """R(g, y) = y - 4 (g - 0.1)^2: tip at (0.1, 0), R_g = 0, R_y = 1, R_gg = -8."""

    def r(g, y):
        return y - 4.0 * (g - 0.1) ** 2

    # lines above the tip cross twice, the line through the tip touches, below: none
    gs = np.linspace(0.0, 0.2, 2001)

    def n_cross(y):
        values = r(gs, y)
        return int(np.sum(np.sign(values[1:]) != np.sign(values[:-1])))

    assert n_cross(0.01) == 2
    assert n_cross(-0.01) == 0
    # at the tip the path derivative vanishes while the in-plane gradient does not
    h = 1e-6
    assert abs((r(0.1 + h, 0.0) - r(0.1 - h, 0.0)) / (2 * h)) < 1e-6
    assert abs((r(0.1, h) - r(0.1, -h)) / (2 * h) - 1.0) < 1e-6
