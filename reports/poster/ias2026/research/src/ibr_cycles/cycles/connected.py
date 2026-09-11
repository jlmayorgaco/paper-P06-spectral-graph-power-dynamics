# ruff: noqa: E501  -- formulas in docstrings and test labels kept on one line
"""Finite-amplitude connected (partition-lattice) cumulants of a characteristic set function.

For a set function ``F`` on the subsets of a finite ground set with ``F(empty) = 1``
(here the collective characteristic factor ``F_s(S) = det(I + Q_SS(s))``), the
connected cumulant of a non-empty ``S`` is

    chi_S = sum_{pi in Pi(S)} (-1)^(|pi|-1) (|pi|-1)!  prod_{B in pi} F(B),

and the moment-cumulant inversion is ``F(S) = sum_{pi in Pi(S)} prod_{B in pi} chi_B``
(theory/FINITE_AMPLITUDE_CONNECTED_CUMULANTS.md, C1).

This is NOT the infinitesimal log-det cumulant of ``cumulants.py``. That object is a
mixed derivative of ``log det(I + Theta K)`` at ``Theta = 0`` (infinitesimal
activation of the unnormalized ``K``). ``chi`` is built from the finite values
``F(B)`` of fully replaced sub-portfolios, is polynomial in them, and stays finite
exactly where ``F`` vanishes (a stability boundary), where any log-based coefficient
is undefined. For scalar ports the two coincide with the multilinear coefficient of
``log det(I + Theta Q)``; for 2x2 ports they differ (theory note, C3).

Keys are ``frozenset`` of labels (bus numbers or indices).
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from itertools import combinations, permutations
from math import factorial

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]
SetFunction = dict[frozenset, complex]


# ----------------------------------------------------------------------- lattices --
def subsets(items: Sequence, *, nonempty: bool = False) -> Iterator[tuple]:
    """All subsets of ``items`` as tuples, by increasing size."""

    start = 1 if nonempty else 0
    for r in range(start, len(items) + 1):
        yield from combinations(items, r)


def set_partitions(items: Sequence) -> Iterator[list[tuple]]:
    """Every set partition of ``items`` exactly once (Bell number many)."""

    items = tuple(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for partition in set_partitions(rest):
        yield [(first,), *partition]
        for i in range(len(partition)):
            yield [*partition[:i], (first, *partition[i]), *partition[i + 1 :]]


def partition_weight(blocks: int) -> int:
    """Moebius function of the partition lattice from a ``blocks``-block partition
    to the one-block partition: ``(-1)^(k-1) (k-1)!``."""

    return (-1) ** (blocks - 1) * factorial(blocks - 1)


# -------------------------------------------------------------- cumulant transforms --
def cumulants_direct(
    values: Mapping[frozenset, complex], ground: Sequence
) -> SetFunction:
    """chi_S by the explicit partition sum (definition). Reference implementation."""

    chi: SetFunction = {}
    for s in subsets(tuple(ground), nonempty=True):
        total = 0.0 + 0.0j
        for pi in set_partitions(s):
            product = 1.0 + 0.0j
            for block in pi:
                product *= values[frozenset(block)]
            total += partition_weight(len(pi)) * product
        chi[frozenset(s)] = complex(total)
    return chi


def cumulants(values: Mapping[frozenset, complex], ground: Sequence) -> SetFunction:
    """chi_S by the recursion ``F(S) = sum_{B: min(S) in B subseteq S} chi_B F(S minus B)``.

    Exact and equivalent to ``cumulants_direct`` (theory note, Lemma 1.2); used for
    speed. ``values`` must contain ``F(empty) = 1``.
    """

    chi: SetFunction = {}
    for s in subsets(tuple(ground), nonempty=True):
        key = frozenset(s)
        anchor = min(s)
        rest = [x for x in s if x != anchor]
        total = values[key]
        for r in range(len(rest)):  # proper blocks containing the anchor
            for extra in combinations(rest, r):
                block = frozenset((anchor, *extra))
                total -= chi[block] * values[key - block]
        chi[key] = complex(total)
    return chi


def moments(chi: Mapping[frozenset, complex], ground: Sequence) -> SetFunction:
    """F(S) = sum_{pi in Pi(S)} prod chi_B, with F(empty) = 1 (C1, forward direction)."""

    out: SetFunction = {frozenset(): 1.0 + 0.0j}
    for s in subsets(tuple(ground), nonempty=True):
        total = 0.0 + 0.0j
        for pi in set_partitions(s):
            product = 1.0 + 0.0j
            for block in pi:
                product *= chi[frozenset(block)]
            total += product
        out[frozenset(s)] = complex(total)
    return out


def boolean_mobius(
    values: Mapping[frozenset, complex], ground: Sequence
) -> SetFunction:
    """mu_T = sum_{R subseteq T} (-1)^(|T|-|R|) F(R): the Boolean interaction terms."""

    out: SetFunction = {}
    for t in subsets(tuple(ground)):
        out[frozenset(t)] = complex(
            sum((-1) ** (len(t) - len(r)) * values[frozenset(r)] for r in subsets(t))
        )
    return out


def partition_terms(
    chi: Mapping[frozenset, complex], s: Sequence
) -> list[tuple[list[tuple], complex]]:
    """Every term ``t_pi = prod_{B in pi} chi_B`` of the expansion of F(S)."""

    terms = []
    for pi in set_partitions(tuple(s)):
        product = 1.0 + 0.0j
        for block in pi:
            product *= chi[frozenset(block)]
        terms.append((pi, complex(product)))
    return terms


def connected_share(chi: Mapping[frozenset, complex], s: Sequence) -> float:
    """nu_S = |chi_S| / sum_pi |t_pi|.

    At a zero of F(S) the composite terms sum to ``-chi_S``, so ``nu_S <= 1/2``;
    ``nu_S = 1/2`` when the composite terms are aligned and ``chi_S`` alone cancels
    them, ``nu_S -> 0`` when the connected term is negligible (theory note, §5).
    """

    terms = partition_terms(chi, s)
    total = sum(abs(t) for _, t in terms)
    return float(abs(chi[frozenset(s)]) / total) if total > 0 else float("nan")


def composite_part(chi: Mapping[frozenset, complex], s: Sequence) -> complex:
    """C_S = F(S) - chi_S = sum over partitions with at least two blocks."""

    return complex(sum(t for pi, t in partition_terms(chi, s) if len(pi) > 1))


def cumulant_support_connected(
    chi: Mapping[frozenset, complex], s: Sequence, tol: float
) -> bool:
    """Is the hypergraph {B subseteq S: |B| >= 2, |chi_B| > tol} connected and spanning S?"""

    s = tuple(s)
    if len(s) == 1:
        return True
    parent = {x: x for x in s}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for b in subsets(s):
        if len(b) >= 2 and abs(chi[frozenset(b)]) > tol:
            root = find(b[0])
            for x in b[1:]:
                parent[find(x)] = root
    return len({find(x) for x in s}) == 1


def algebraic_order(mu: Mapping[frozenset, complex], tol: float) -> int:
    """d_alg = max{|T| : |mu_T| > tol} (0 if only the empty set survives)."""

    return max((len(t) for t, v in mu.items() if abs(v) > tol), default=0)


def connected_order(chi: Mapping[frozenset, complex], tol: float) -> int:
    """d_conn = max{|T| >= 2 : |chi_T| > tol}, or 1 if no multi-element cumulant survives."""

    return max(
        (len(t) for t, v in chi.items() if len(t) >= 2 and abs(v) > tol), default=1
    )


# --------------------------------------------------------- characteristic functions --
def channel_index(block_sizes: Sequence[int]) -> list[list[int]]:
    out, start = [], 0
    for size in block_sizes:
        out.append(list(range(start, start + size)))
        start += size
    return out


def characteristic_values(q: ComplexMatrix, block_sizes: Sequence[int]) -> SetFunction:
    """F(S) = det(I + Q_SS) over block index sets, F(empty) = 1. Keys are block indices."""

    chans = channel_index(block_sizes)
    out: SetFunction = {}
    for s in subsets(tuple(range(len(block_sizes)))):
        cols = [c for b in s for c in chans[b]]
        out[frozenset(s)] = (
            complex(np.linalg.det(np.eye(len(cols)) + q[np.ix_(cols, cols)]))
            if cols
            else 1.0 + 0.0j
        )
    return out


def spanning_cycle_sum(q: ComplexMatrix, s: Sequence[int]) -> complex:
    """Sum over the (|S|-1)! spanning directed cycles gamma of S of prod_{(i,j) in gamma} q_ij."""

    s = tuple(s)
    head, rest = s[0], s[1:]
    total = 0.0 + 0.0j
    for tail in permutations(rest):
        cycle = (head, *tail)
        product = 1.0 + 0.0j
        for k, i in enumerate(cycle):
            product *= q[i, cycle[(k + 1) % len(cycle)]]
        total += product
    return complex(total)


def _cycles_of(perm: dict[int, int]) -> list[list[int]]:
    seen, cycles = set(), []
    for start in perm:
        if start in seen:
            continue
        cycle, x = [], start
        while x not in seen:
            seen.add(x)
            cycle.append(x)
            x = perm[x]
        cycles.append(cycle)
    return cycles


def block_connected_expansion(
    q: ComplexMatrix, block_sizes: Sequence[int], s: Sequence[int]
) -> complex:
    """Brute-force sum over the block-connected permutations of the channels of S.

    sigma is block-connected when the blocks of S are connected by the relation
    "some cycle of sigma contains channels of both". Each sigma contributes
    ``sgn(sigma) prod_c (I + Q)_{c, sigma(c)}`` (theory note, Theorem 3.3).
    Exponential cost; for verification only (|channels| <= 8).
    """

    chans = channel_index(block_sizes)
    owner = {c: b for b in s for c in chans[b]}
    cols = [c for b in s for c in chans[b]]
    ipq = np.eye(q.shape[0], dtype=np.complex128) + q
    total = 0.0 + 0.0j
    for image in permutations(cols):
        perm = dict(zip(cols, image, strict=True))
        cycles = _cycles_of(perm)
        parent = {b: b for b in s}

        def find(x, parent=parent):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for cycle in cycles:
            blocks = {owner[c] for c in cycle}
            first = next(iter(blocks))
            for b in blocks:
                parent[find(b)] = find(first)
        if len({find(b) for b in s}) != 1:
            continue
        sign = 1
        for cycle in cycles:
            sign *= (-1) ** (len(cycle) - 1)
        product = 1.0 + 0.0j
        for c in cols:
            product *= ipq[c, perm[c]]
        total += sign * product
    return complex(total)


def block_cycle_trace_sum(
    q: ComplexMatrix, block_sizes: Sequence[int], s: Sequence[int]
) -> complex:
    """Sum over cyclic orders b1 -> ... -> bk -> b1 of S of tr(Q_b1b2 ... Q_bkb1).

    The single-cycle, one-visit-per-block part of the block expansion; equal to the
    multilinear coefficient of log det(I + Theta Q) with Theta = blkdiag(theta_b I)
    (the naive 'holonomy' formula). Not equal to chi_S for 2x2 ports.
    """

    chans = channel_index(block_sizes)
    s = tuple(s)
    head, rest = s[0], s[1:]
    total = 0.0 + 0.0j
    for tail in permutations(rest):
        cycle = (head, *tail)
        product = np.eye(len(chans[cycle[0]]), dtype=np.complex128)
        for k, b in enumerate(cycle):
            nxt = cycle[(k + 1) % len(cycle)]
            product = product @ q[np.ix_(chans[b], chans[nxt])]
        total += np.trace(product)
    return complex(total)


# ------------------------------------------------------------------- derivatives --
def adjugate(x: ComplexMatrix) -> ComplexMatrix:
    """Classical adjugate by cofactors (valid for singular x; small matrices)."""

    n = x.shape[0]
    if n == 1:
        return np.ones((1, 1), dtype=np.complex128)
    adj = np.empty((n, n), dtype=np.complex128)
    for i in range(n):
        for j in range(n):
            minor = np.delete(np.delete(x, i, axis=0), j, axis=1)
            adj[j, i] = (-1) ** (i + j) * np.linalg.det(minor)
    return adj


def characteristic_derivatives(
    q: ComplexMatrix, dq: ComplexMatrix, block_sizes: Sequence[int]
) -> SetFunction:
    """dF(S) = tr(adj(I + Q_SS) dQ_SS) (Jacobi's formula in adjugate form; exact at a zero)."""

    chans = channel_index(block_sizes)
    out: SetFunction = {}
    for s in subsets(tuple(range(len(block_sizes)))):
        cols = [c for b in s for c in chans[b]]
        if not cols:
            out[frozenset(s)] = 0.0 + 0.0j
            continue
        sub = np.eye(len(cols)) + q[np.ix_(cols, cols)]
        out[frozenset(s)] = complex(np.trace(adjugate(sub) @ dq[np.ix_(cols, cols)]))
    return out


def cumulant_derivatives(
    values: Mapping[frozenset, complex],
    dvalues: Mapping[frozenset, complex],
    ground: Sequence,
) -> SetFunction:
    """d chi_S from dF by differentiating the partition formula (product rule over blocks)."""

    out: SetFunction = {}
    for s in subsets(tuple(ground), nonempty=True):
        total = 0.0 + 0.0j
        for pi in set_partitions(s):
            vals = [values[frozenset(b)] for b in pi]
            dvals = [dvalues[frozenset(b)] for b in pi]
            for k in range(len(pi)):
                product = dvals[k]
                for j in range(len(pi)):
                    if j != k:
                        product *= vals[j]
                total += partition_weight(len(pi)) * product
        out[frozenset(s)] = complex(total)
    return out
