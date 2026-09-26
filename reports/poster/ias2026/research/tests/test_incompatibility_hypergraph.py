"""F2C: the minimal spectral incompatibility hypergraph.

Formal content of ``theory/F2C_incompatibility_hypergraph.md``. The random tests
check the combinatorial statements on arbitrary count assignments; the parametric
test is the explicit construction behind "same kappa, different hypergraph".
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pytest

from ibr_cycles.diagnosis.composability import (
    NO_FINITE_ORDER,
    BandRightHalfPlane,
    composability_order,
    hypergraph_label,
    hypergraph_move,
    hypergraph_order,
    incompatibility_hypergraph,
    is_antichain,
    maximal_free_sets,
    parse_hypergraph_label,
    spectral_count,
    upward_closure_defect,
)

BAND = BandRightHalfPlane(low_hz=0.3, high_hz=1.5)


def power_set(ground):
    return [
        frozenset(c) for k in range(len(ground) + 1) for c in combinations(ground, k)
    ]


def random_counts(rng, ground, p_unsafe=0.35):
    counts = {
        s: int(rng.random() < p_unsafe) * int(rng.integers(1, 3))
        for s in power_set(ground)
    }
    counts[frozenset()] = 0
    return counts


@pytest.mark.parametrize("seed", range(200))
def test_hypergraph_is_an_antichain_and_carries_kappa(seed):
    rng = np.random.default_rng(seed)
    counts = random_counts(rng, (1, 2, 3, 4))
    edges = incompatibility_hypergraph(counts)
    assert is_antichain(edges)
    assert hypergraph_order(edges) == composability_order(counts)
    # every unsafe set contains a hyperedge, and every hyperedge is unsafe
    for s, n in counts.items():
        if s and n > 0:
            assert any(e <= s for e in edges)
    for e in edges:
        assert counts[e] > 0
        assert all(counts[r] == 0 for r in counts if r < e)


@pytest.mark.parametrize("seed", range(50))
def test_hypergraph_restricts_consistently_to_a_sub_fleet(seed):
    """H of the power set of B equals the hyperedges of H(A) contained in B."""

    rng = np.random.default_rng(1000 + seed)
    ground = (1, 2, 3, 4, 5)
    counts = random_counts(rng, ground, p_unsafe=0.25)
    full = set(incompatibility_hypergraph(counts))
    sub = frozenset((1, 2, 4))
    restricted = {s: n for s, n in counts.items() if s <= sub}
    assert set(incompatibility_hypergraph(restricted)) == {e for e in full if e <= sub}


def test_empty_hypergraph_is_infinite_order():
    counts = {s: 0 for s in power_set((1, 2, 3))}
    assert incompatibility_hypergraph(counts) == ()
    assert hypergraph_order(()) == NO_FINITE_ORDER
    assert hypergraph_label(()) == "EMPTY"


def test_label_round_trip_and_canonical_order():
    edges = (frozenset({37, 30}), frozenset({30, 33, 35}), frozenset({33}))
    label = hypergraph_label(edges)
    assert label == "33|30+37|30+33+35"
    assert set(parse_hypergraph_label(label)) == set(edges)


def test_non_downward_closed_family_is_refused():
    counts = {frozenset(): 0, frozenset({1, 2}): 1, frozenset({1}): 0}
    with pytest.raises(ValueError, match="downward closed"):
        incompatibility_hypergraph(counts)


def test_upward_closure_defect_detects_a_restabilised_superset():
    """Dynamic stability is not monotone: a superset of a hyperedge can be safe."""

    counts = {s: 0 for s in power_set((1, 2, 3))}
    counts[frozenset({1, 2})] = 1
    counts[frozenset({1, 2, 3})] = 0
    assert incompatibility_hypergraph(counts) == (frozenset({1, 2}),)
    assert upward_closure_defect(counts) == (frozenset({1, 2, 3}),)


@pytest.mark.parametrize("seed", range(100))
def test_hyperedge_free_portfolios_are_safe_without_heredity(seed):
    rng = np.random.default_rng(5000 + seed)
    ground = (1, 2, 3, 4)
    counts = random_counts(rng, ground)
    edges = incompatibility_hypergraph(counts)
    for s in maximal_free_sets(edges, ground):
        assert all(counts[r] == 0 for r in counts if r <= s)
    # and they are maximal: adding any vertex creates a hyperedge
    for s in maximal_free_sets(edges, ground):
        for v in set(ground) - s:
            assert any(e <= s | {v} for e in edges)


def test_maximal_free_sets_of_the_flagship_hypergraphs():
    core = (30, 33, 35, 37)
    assert maximal_free_sets([frozenset(core)], core) == tuple(
        frozenset(c) for c in [(30, 33, 35), (30, 33, 37), (30, 35, 37), (33, 35, 37)]
    )
    assert set(maximal_free_sets([frozenset({30, 33})], core)) == {
        frozenset({30, 35, 37}),
        frozenset({33, 35, 37}),
    }


def test_contraction_move_is_recognised():
    move = hypergraph_move([frozenset({30, 33, 35, 37})], [frozenset({30, 33, 35})])
    assert move["kind"] == "CONTRACTION"
    move = hypergraph_move([frozenset({30, 33, 35})], [frozenset({30, 33})])
    assert move["kind"] == "CONTRACTION"
    move = hypergraph_move([frozenset({1, 2})], [frozenset({1, 2}), frozenset({1, 3})])
    assert move["kind"] == "INSERTION"


# ------------------------------------------ same kappa, different hypergraph --

OMEGA0 = 2.0 * np.pi * 0.6


def toy_matrix(members, theta):
    """One oscillatory mode per replacement set, real part additive in the actions.

    sigma(S; theta) = -1 + 0.6 [1 in S] + (0.6 - 0.4 theta) [2 in S]
                         + (0.3 + 0.4 theta) [3 in S]
    """

    sigma = -1.0
    sigma += 0.6 * (1 in members)
    sigma += (0.6 - 0.4 * theta) * (2 in members)
    sigma += (0.3 + 0.4 * theta) * (3 in members)
    return np.array([[sigma, OMEGA0], [-OMEGA0, sigma]])


def toy_hypergraph(theta):
    counts = {
        s: spectral_count(np.linalg.eigvals(toy_matrix(s, theta)), BAND)
        for s in power_set((1, 2, 3))
    }
    return composability_order(counts), incompatibility_hypergraph(counts)


def test_same_kappa_different_hypergraph_on_a_continuous_path():
    """kappa stays 2 while H goes {12} -> {12, 13} -> {13}; boundaries 0.25, 0.5."""

    seen = {}
    for theta in np.linspace(0.0, 1.0, 41):
        if min(abs(theta - 0.25), abs(theta - 0.5)) < 1e-9:
            continue  # exactly on a spectral boundary: outside Theta_reg
        kappa, edges = toy_hypergraph(theta)
        assert kappa == 2
        seen.setdefault(hypergraph_label(edges), []).append(theta)
    assert set(seen) == {"1+2", "1+2|1+3", "1+3"}
    assert max(seen["1+2"]) < 0.25 < min(seen["1+2|1+3"])
    assert max(seen["1+2|1+3"]) < 0.5 < min(seen["1+3"])
