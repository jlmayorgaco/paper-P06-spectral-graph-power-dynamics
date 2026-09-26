# ruff: noqa: E501  -- formulas in docstrings and test labels kept on one line
"""Finite-amplitude connected cumulants: C1-C4 of the theory note and toy tests T1-T5.

theory/FINITE_AMPLITUDE_CONNECTED_CUMULANTS.md. Random inputs use fixed seeds; no
example is tuned except where a test states the tuning (C4 (d)).
"""

from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np
import pytest

from ibr_cycles.cycles.connected import (
    algebraic_order,
    block_connected_expansion,
    block_cycle_trace_sum,
    boolean_mobius,
    channel_index,
    characteristic_derivatives,
    characteristic_values,
    composite_part,
    connected_order,
    connected_share,
    cumulant_derivatives,
    cumulant_support_connected,
    cumulants,
    cumulants_direct,
    moments,
    set_partitions,
    spanning_cycle_sum,
    subsets,
)
from ibr_cycles.models.lag_network_toy import (
    LagNetworkToy,
    path,
    ring,
    ring_plus_pair,
    scaled_to_boundary,
    star,
    tuned_mobius_zero,
    two_pairs,
)

BELL = [1, 1, 2, 5, 15, 52, 203, 877]
TOL = 1e-11


def rand_c(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def random_setfunction(rng, n):
    f = {frozenset(s): complex(rand_c(rng)) for s in subsets(range(n), nonempty=True)}
    f[frozenset()] = 1.0 + 0.0j
    return f


def zero_diag_blocks(q, p):
    q = q.copy()
    for b in range(q.shape[0] // p):
        q[p * b : p * b + p, p * b : p * b + p] = 0.0
    return q


# ------------------------------------------------------------- lattice machinery --
@pytest.mark.parametrize("n", range(8))
def test_bell_numbers(n):
    parts = list(set_partitions(tuple(range(n))))
    assert len(parts) == BELL[n]
    assert len({tuple(sorted(tuple(sorted(b)) for b in p)) for p in parts}) == BELL[n]


@pytest.mark.parametrize("n", range(1, 7))
def test_recursion_equals_partition_definition(n):
    rng = np.random.default_rng(100 + n)
    f = random_setfunction(rng, n)
    a, b = cumulants(f, range(n)), cumulants_direct(f, range(n))
    assert max(abs(a[k] - b[k]) for k in a) < 1e-10 * max(
        1, max(abs(v) for v in b.values())
    )


# ---------------------------------------------------------------------------- C1 --
@pytest.mark.parametrize("n", range(1, 7))
def test_C1_moment_cumulant_inversion(n):
    rng = np.random.default_rng(200 + n)
    f = random_setfunction(rng, n)
    back = moments(cumulants(f, range(n)), range(n))
    assert max(abs(back[k] - f[k]) for k in f) < 1e-10 * max(abs(v) for v in f.values())


@pytest.mark.parametrize("n", range(1, 7))
def test_C1_uniqueness_triangularity(n):
    """Cumulants -> moments -> cumulants is the identity (inverse on both sides)."""

    rng = np.random.default_rng(250 + n)
    chi = {frozenset(s): complex(rand_c(rng)) for s in subsets(range(n), nonempty=True)}
    again = cumulants(moments(chi, range(n)), range(n))
    assert max(abs(again[k] - chi[k]) for k in chi) < 1e-9


# ---------------------------------------------------------------------------- C2 --
@pytest.mark.parametrize("split", [(2, 2), (1, 3), (3, 3), (2, 4)])
def test_C2_factorization_annihilates_every_mixed_cumulant(split):
    n1, n2 = split
    rng = np.random.default_rng(300 + 10 * n1 + n2)
    g1, g2 = random_setfunction(rng, n1), random_setfunction(rng, n2)
    s1, s2 = set(range(n1)), set(range(n1, n1 + n2))
    f = {}
    for a in subsets(range(n1 + n2)):
        a1 = frozenset(x for x in a if x in s1)
        a2 = frozenset(x - n1 for x in a if x in s2)
        f[frozenset(a)] = g1[a1] * g2[a2]
    chi = cumulants(f, range(n1 + n2))
    mixed = [abs(v) for k, v in chi.items() if k & s1 and k & s2]
    assert max(mixed) < 1e-10
    # and the unmixed cumulants are those of the factors
    for k, v in chi.items():
        if k <= s1:
            assert abs(v - cumulants(g1, range(n1))[k]) < 1e-10


def test_C2_converse_vanishing_mixed_cumulants_gives_factorization():
    rng = np.random.default_rng(310)
    n = 5
    s1 = {0, 1}
    chi = {}
    for k in subsets(range(n), nonempty=True):
        key = frozenset(k)
        chi[key] = 0j if (key & s1 and key - s1) else complex(rand_c(rng))
    f = moments(chi, range(n))
    for a in subsets(range(n)):
        a = frozenset(a)
        assert abs(f[a] - f[a & s1] * f[a - s1]) < 1e-10


def test_C2_top_cumulant_zero_does_not_imply_factorization():
    """Path coupling 0-1-2-3: chi_0123 = 0 but F(0123) != F(01) F(23) etc."""

    q = np.zeros((4, 4), complex)
    for i in range(3):
        q[i, i + 1], q[i + 1, i] = 0.7, 0.9
    f = characteristic_values(q, [1] * 4)
    chi = cumulants(f, range(4))
    assert abs(chi[frozenset(range(4))]) < TOL
    full = f[frozenset(range(4))]
    for a in [(0,), (0, 1), (0, 1, 2)]:
        a = frozenset(a)
        rest = frozenset(range(4)) - a
        assert abs(full - f[a] * f[rest]) > 1e-3


# ---------------------------------------------------------------------------- C3 --
@pytest.mark.parametrize("n", range(2, 7))
@pytest.mark.parametrize("zero_diagonal", [True, False])
def test_C3_scalar_spanning_cycle_formula(n, zero_diagonal):
    rng = np.random.default_rng(400 + n + 50 * zero_diagonal)
    q = rand_c(rng, n, n) * 0.6
    if zero_diagonal:
        np.fill_diagonal(q, 0.0)
    chi = cumulants(characteristic_values(q, [1] * n), range(n))
    for s in subsets(range(n), nonempty=True):
        if len(s) == 1:
            assert abs(chi[frozenset(s)] - (1.0 + q[s[0], s[0]])) < TOL
        else:
            expected = (-1) ** (len(s) - 1) * spanning_cycle_sum(q, s)
            assert abs(chi[frozenset(s)] - expected) < 1e-10 * max(1.0, abs(expected))


def _lambda_log(poly, ground):
    """log of 1 + nilpotent in the algebra C[x]/(x_i^2): dict frozenset -> coefficient."""

    def mul(p1, p2):
        out = {}
        for k1, v1 in p1.items():
            for k2, v2 in p2.items():
                if k1 & k2:
                    continue
                out[k1 | k2] = out.get(k1 | k2, 0) + v1 * v2
        return out

    u = {k: v for k, v in poly.items() if k}
    result, power = {}, {frozenset(): 1.0 + 0.0j}
    for k in range(1, len(ground) + 1):
        power = mul(power, u)
        for key, v in power.items():
            result[key] = result.get(key, 0) + (-1) ** (k - 1) / k * v
    return result


@pytest.mark.parametrize("n", range(2, 6))
def test_C3_scalar_cumulant_equals_multilinear_logdet_coefficient(n):
    """Scalar ports: chi_S = [theta^S] log det(I + diag(theta) Q)."""

    rng = np.random.default_rng(450 + n)
    q = rand_c(rng, n, n) * 0.5
    np.fill_diagonal(q, 0.0)
    poly = {
        frozenset(u): complex(np.linalg.det(q[np.ix_(u, u)])) if u else 1.0 + 0j
        for u in subsets(range(n))
    }
    log_poly = _lambda_log(poly, range(n))
    chi = cumulants(characteristic_values(q, [1] * n), range(n))
    for s in subsets(range(n)):
        if len(s) >= 2:
            assert abs(log_poly.get(frozenset(s), 0) - chi[frozenset(s)]) < 1e-10


@pytest.mark.parametrize("m", [2, 3])
def test_C3_block_connected_expansion_is_exact(m):
    rng = np.random.default_rng(500 + m)
    q = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.5, 2)
    chi = cumulants(characteristic_values(q, [2] * m), range(m))
    for s in subsets(range(m), nonempty=True):
        assert abs(chi[frozenset(s)] - block_connected_expansion(q, [2] * m, s)) < 1e-10


def test_C3_block_connected_expansion_four_blocks():
    rng = np.random.default_rng(504)
    q = zero_diag_blocks(rand_c(rng, 8, 8) * 0.5, 2)
    chi = cumulants(characteristic_values(q, [2] * 4), range(4))
    assert (
        abs(chi[frozenset(range(4))] - block_connected_expansion(q, [2] * 4, range(4)))
        < 1e-10
    )


def test_C3_two_block_closed_form_and_naive_formula_refuted():
    rng = np.random.default_rng(510)
    q = zero_diag_blocks(rand_c(rng, 4, 4) * 0.8, 2)
    chi = cumulants(characteristic_values(q, [2, 2]), range(2))[frozenset({0, 1})]
    x, y = q[0:2, 2:4], q[2:4, 0:2]
    exact = -np.trace(x @ y) + np.linalg.det(x) * np.linalg.det(y)
    naive = -np.trace(x @ y)  # the scalar cycle formula with q_ij -> holonomy trace
    assert abs(chi - exact) < 1e-12
    assert abs(chi - naive) > 1e-2


@pytest.mark.parametrize("m", [2, 3, 4])
def test_C3_block_single_cycle_part_equals_multilinear_logdet(m):
    """For 2x2 ports the naive holonomy sum is the infinitesimal cumulant, not chi."""

    rng = np.random.default_rng(520 + m)
    q = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.5, 2)
    chans = channel_index([2] * m)
    poly = {}
    for u in subsets(range(2 * m)):
        blocks = [c // 2 for c in u]
        if len(set(blocks)) != len(blocks):
            continue  # theta_b^2 = 0
        poly[frozenset(blocks)] = poly.get(frozenset(blocks), 0) + (
            complex(np.linalg.det(q[np.ix_(u, u)])) if u else 1.0 + 0j
        )
    log_poly = _lambda_log(poly, range(m))
    chi = cumulants(characteristic_values(q, [2] * m), range(m))
    assert len(chans) == m
    for s in subsets(range(m)):
        if len(s) >= 2:
            naive = (-1) ** (len(s) - 1) * block_cycle_trace_sum(q, [2] * m, s)
            assert abs(log_poly.get(frozenset(s), 0) - naive) < 1e-10
            assert abs(chi[frozenset(s)] - naive) > 1e-6


def test_block_gauge_invariance():
    rng = np.random.default_rng(530)
    m = 3
    q = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.5, 2)
    s_b = np.zeros((2 * m, 2 * m), complex)
    for b in range(m):
        s_b[2 * b : 2 * b + 2, 2 * b : 2 * b + 2] = rand_c(rng, 2, 2)
    q2 = np.linalg.solve(s_b, q @ s_b)
    c1 = cumulants(characteristic_values(q, [2] * m), range(m))
    c2 = cumulants(characteristic_values(q2, [2] * m), range(m))
    assert max(abs(c1[k] - c2[k]) for k in c1) < 1e-10
    for s in [(0, 1), (0, 1, 2)]:
        assert (
            abs(
                block_cycle_trace_sum(q, [2] * m, s)
                - block_cycle_trace_sum(q2, [2] * m, s)
            )
            < 1e-10
        )


# --------------------------------------------------------- Moebius <-> cumulants --
@pytest.mark.parametrize("n", range(1, 6))
def test_mobius_is_partition_sum_of_reduced_cumulants(n):
    """mu_S = sum_pi prod psi_B, psi = chi except psi_{i} = chi_{i} - 1."""

    rng = np.random.default_rng(600 + n)
    f = random_setfunction(rng, n)
    chi = cumulants(f, range(n))
    psi = {k: (v - 1.0 if len(k) == 1 else v) for k, v in chi.items()}
    mu = boolean_mobius(f, range(n))
    for s in subsets(range(n), nonempty=True):
        total = 0j
        for pi in set_partitions(s):
            t = 1 + 0j
            for b in pi:
                t *= psi[frozenset(b)]
            total += t
        assert abs(total - mu[frozenset(s)]) < 1e-9


# ------------------------------------------------------------------- derivatives --
def test_derivatives_match_finite_differences_and_hold_at_a_zero():
    rng = np.random.default_rng(700)
    m = 3
    q0 = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.5, 2)
    q1 = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.5, 2)
    # rescale so that Q has an eigenvalue at -1: det(I + Q_VV) = 0 (a characteristic zero)
    ev = np.linalg.eigvals(q0)
    q0 = q0 / (-ev[np.argmax(np.abs(ev))])
    h = 1e-6
    f0 = characteristic_values(q0, [2] * m)
    assert abs(f0[frozenset(range(m))]) < 1e-10
    df = characteristic_derivatives(q0, q1, [2] * m)
    fp = characteristic_values(q0 + h * q1, [2] * m)
    fm = characteristic_values(q0 - h * q1, [2] * m)
    for k in f0:
        assert abs(df[k] - (fp[k] - fm[k]) / (2 * h)) < 1e-6
    dchi = cumulant_derivatives(f0, df, range(m))
    cp, cm = cumulants(fp, range(m)), cumulants(fm, range(m))
    for k in dchi:
        assert abs(dchi[k] - (cp[k] - cm[k]) / (2 * h)) < 1e-6


# --------------------------------------------------------------- toys T1 - T5 --
S_PROBES = [0.1 + 0.4j, 0.0 + 1.3j, 0.5 - 0.2j, 2.0 + 0.0j]


def test_T1_independent_pairs_full_mobius_nonzero_but_top_cumulant_zero():
    toy = two_pairs()
    for s in S_PROBES:
        f = characteristic_values(toy.q(s), [1] * 4)
        mu = boolean_mobius(f, range(4))
        chi = cumulants(f, range(4))
        full = frozenset(range(4))
        assert abs(mu[full]) > 1e-3
        assert abs(chi[full]) < TOL
        assert abs(mu[full] - chi[frozenset({0, 1})] * chi[frozenset({2, 3})]) < TOL


def test_T2_ring_top_cumulant_is_the_spanning_cycle():
    toy = ring()
    for s in S_PROBES:
        q = toy.q(s)
        chi = cumulants(characteristic_values(q, [1] * 4), range(4))
        expected = -q[0, 1] * q[1, 2] * q[2, 3] * q[3, 0]
        assert abs(chi[frozenset(range(4))] - expected) < TOL
        assert abs(chi[frozenset(range(4))] + spanning_cycle_sum(q, range(4))) < TOL
        for b in subsets(range(4)):
            if 2 <= len(b) < 4:
                assert abs(chi[frozenset(b)]) < TOL


def boundary_profile(toy):
    target = tuple(range(toy.m))
    at = scaled_to_boundary(toy, target)
    s_star = complex(0.0, abs(at.critical(target).imag))
    # zero of F(V) at the boundary eigenvalue
    lam = at.critical(target)
    f = characteristic_values(at.q(lam), [toy.block] * toy.m)
    chi = cumulants(f, range(toy.m))
    return at, lam, s_star, f, chi


def test_T3_same_kappa_four_composite_versus_connected():
    # composite: symmetric path (chi_V identically 0); connected: directed ring
    for maker, expect_nu in ((path, 0.0), (ring, 0.5)):
        at, lam, _, f, chi = boundary_profile(maker())
        full = frozenset(range(4))
        assert abs(f[full]) < 1e-8  # at the characteristic zero
        # proper subsets strictly stable at the boundary, full set marginal
        assert abs(at.alpha(tuple(range(4)))) < 1e-9
        for r in (1, 2, 3):
            for s in combinations(range(4), r):
                assert at.alpha(s) < -1e-6
        assert connected_share(chi, range(4)) == pytest.approx(expect_nu, abs=1e-8)
        assert cumulant_support_connected(chi, range(4), 1e-10)
        assert abs(chi[full] + composite_part(chi, range(4))) < 1e-8


def test_T3_kappa_four_is_realized_away_from_the_boundary():
    for maker in (path, ring, star):
        at = scaled_to_boundary(maker(), (0, 1, 2, 3))
        pushed = LagNetworkToy(1.02 * at.g, at.a)
        assert pushed.hypergraph() == [(0, 1, 2, 3)]
        assert pushed.kappa() == 4


@pytest.mark.parametrize("n", range(2, 7))
def test_T4_random_scalar_partition_vs_cycle_enumeration(n):
    rng = np.random.default_rng(800 + n)
    for _ in range(5):
        q = rand_c(rng, n, n) * 0.7
        np.fill_diagonal(q, 0.0)
        chi = cumulants_direct(characteristic_values(q, [1] * n), range(n))
        for s in subsets(range(n)):
            if len(s) >= 2:
                assert (
                    abs(
                        chi[frozenset(s)]
                        - (-1) ** (len(s) - 1) * spanning_cycle_sum(q, s)
                    )
                    < 1e-10
                )


@pytest.mark.parametrize("m", [2, 3, 4])
def test_T5_random_blocks_inversion_and_factorization(m):
    rng = np.random.default_rng(900 + m)
    q = zero_diag_blocks(rand_c(rng, 2 * m, 2 * m) * 0.6, 2)
    f = characteristic_values(q, [2] * m)
    chi = cumulants(f, range(m))
    back = moments(chi, range(m))
    assert max(abs(back[k] - f[k]) for k in f) < 1e-10
    # decouple blocks {0} and {1..m-1}: zero the cross blocks -> factorization -> mixed chi = 0
    qd = q.copy()
    qd[0:2, 2:] = 0.0
    qd[2:, 0:2] = 0.0
    chi_d = cumulants(characteristic_values(qd, [2] * m), range(m))
    assert max(abs(v) for k, v in chi_d.items() if 0 in k and len(k) > 1) < 1e-12


# ---------------------------------------------------------------------------- C4 --
def orders(toy, tol=1e-9):
    d_alg = d_conn = 0
    for s in S_PROBES:
        f = characteristic_values(toy.q(s), [toy.block] * toy.m)
        d_alg = max(d_alg, algebraic_order(boolean_mobius(f, range(toy.m)), tol))
        d_conn = max(d_conn, connected_order(cumulants(f, range(toy.m)), tol))
    return d_alg, d_conn


def test_C4_counterexamples():
    ring_b = LagNetworkToy(1.02 * scaled_to_boundary(ring(), (0, 1, 2, 3)).g, 1.0)
    path_b = LagNetworkToy(1.02 * scaled_to_boundary(path(), (0, 1, 2, 3)).g, 1.0)
    star_b = LagNetworkToy(1.02 * scaled_to_boundary(star(), (0, 1, 2, 3)).g, 1.0)
    assert (ring_b.kappa(), *orders(ring_b)) == (4, 4, 4)
    assert (path_b.kappa(), *orders(path_b)) == (
        4,
        4,
        2,
    )  # kappa > d_conn, d_alg > d_conn
    assert (star_b.kappa(), *orders(star_b)) == (4, 2, 2)  # kappa > d_alg
    tuned = tuned_mobius_zero()
    d_alg, d_conn = orders(tuned)
    assert d_conn == 4 and d_alg < 4  # d_conn > d_alg
    rp = ring_plus_pair()
    assert rp.kappa() == 2
    assert orders(rp) == (4, 4)  # kappa < d_alg = d_conn
    assert comb(4, 2) == 6


def test_minimal_coalitions_are_cumulant_connected_in_random_toys():
    """Theorem 2.3: at the characteristic zero of a minimal coalition the cumulant
    support is connected and spans it. Random toys, no tuning."""

    rng = np.random.default_rng(1000)
    checked = 0
    for _ in range(200):
        g = rng.normal(size=(4, 4)) * (rng.random((4, 4)) < 0.5)
        np.fill_diagonal(g, 0.0)
        toy = LagNetworkToy(g.astype(complex), 1.0)
        for h in toy.hypergraph():
            if len(h) < 2:
                continue
            lam = toy.critical(h)
            sub = LagNetworkToy(g[np.ix_(h, h)].astype(complex), 1.0)
            f = characteristic_values(sub.q(lam), [1] * len(h))
            chi = cumulants(f, range(len(h)))
            assert abs(f[frozenset(range(len(h)))]) < 1e-8
            assert cumulant_support_connected(chi, range(len(h)), 1e-9)
            checked += 1
    assert checked > 20
