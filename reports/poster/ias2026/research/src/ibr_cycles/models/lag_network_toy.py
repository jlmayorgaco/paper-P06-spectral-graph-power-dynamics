# ruff: noqa: E501  -- formulas in docstrings and test labels kept on one line
"""Lag-network toy family for the connected-cumulant counterexamples (synthetic).

Ports i = 1..m, each a block of ``p`` channels. Replacing the set S switches on the
units of S, whose dynamics are

    x_S' = A_S x_S,    A_S = -a I - G_SS,

with G block-zero-diagonal. The normalized port interaction is

    Q(s) = G / (s + a),     F_s(S) = det(I + Q_SS(s)) = det(sI - A_S) / (s + a)^(p|S|).

For a > 0 the base and every singleton are stable (A_i = -a I), F_s(S) is analytic
on the closed right half-plane, and its zeros there are exactly the unstable
eigenvalues of A_S: the admissibility hypothesis (A) of the theory note holds
exactly. It is a purely synthetic device, not a power-system model.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class LagNetworkToy:
    g: NDArray[np.complex128]
    a: float
    block: int = 1

    @property
    def m(self) -> int:
        return self.g.shape[0] // self.block

    def cols(self, s) -> list[int]:
        return [self.block * i + c for i in sorted(s) for c in range(self.block)]

    def q(self, s: complex) -> NDArray[np.complex128]:
        return np.asarray(self.g, dtype=np.complex128) / (s + self.a)

    def a_matrix(self, s) -> NDArray[np.complex128]:
        c = self.cols(s)
        return (
            -self.a * np.eye(len(c))
            - np.asarray(self.g, dtype=np.complex128)[np.ix_(c, c)]
        )

    def spectrum(self, s) -> NDArray[np.complex128]:
        return (
            np.linalg.eigvals(self.a_matrix(s))
            if s
            else np.zeros(0, dtype=np.complex128)
        )

    def alpha(self, s) -> float:
        ev = self.spectrum(s)
        return float(ev.real.max()) if ev.size else -np.inf

    def unstable(self, s) -> bool:
        return self.alpha(s) > 0.0

    def hypergraph(self) -> list[tuple[int, ...]]:
        """Minimal unstable portfolios (antichain)."""

        found: list[tuple[int, ...]] = []
        for r in range(1, self.m + 1):
            for s in combinations(range(self.m), r):
                if any(set(h) <= set(s) for h in found):
                    continue
                if self.unstable(s):
                    found.append(s)
        return found

    def kappa(self) -> float:
        h = self.hypergraph()
        return float(min(len(e) for e in h)) if h else float("inf")

    def critical(self, s) -> complex:
        ev = self.spectrum(s)
        return complex(ev[np.argmax(ev.real)])


def two_pairs(a: float = 1.0) -> LagNetworkToy:
    """Two independent reciprocal pairs {0,1} and {2,3}; no cross-coupling (T1)."""

    g = np.zeros((4, 4), complex)
    g[0, 1] = g[1, 0] = g[2, 3] = g[3, 2] = -0.8
    return LagNetworkToy(g, a)


def ring(w: float = 1.0, a: float = 1.0) -> LagNetworkToy:
    """Directed ring 0 -> 1 -> 2 -> 3 -> 0: a single spanning cycle (T2, T3 connected)."""

    g = np.zeros((4, 4), complex)
    for i in range(4):
        g[i, (i + 1) % 4] = w
    return LagNetworkToy(g, a)


def path(w: float = -1.0, a: float = 1.0) -> LagNetworkToy:
    """Symmetric path 0 - 1 - 2 - 3: only pair clusters (T3 composite)."""

    g = np.zeros((4, 4), complex)
    for i in range(3):
        g[i, i + 1] = g[i + 1, i] = w
    return LagNetworkToy(g, a)


def star(w: float = -1.0, a: float = 1.0) -> LagNetworkToy:
    """Symmetric star with centre 0 (singular adjacency: mu_V identically 0)."""

    g = np.zeros((4, 4), complex)
    for leaf in (1, 2, 3):
        g[0, leaf] = g[leaf, 0] = w
    return LagNetworkToy(g, a)


def tuned_mobius_zero() -> LagNetworkToy:
    """Pairs 0-1, 2-3 plus the ring 0 -> 2 -> 1 -> 3 -> 0, with the ring weight
    g[3, 0] solved (det is affine in it) so that det G = 0: mu_V identically 0
    while chi_V != 0. Tuned by construction, used only as a counterexample (C4 d)."""

    g = np.zeros((4, 4), complex)
    g[0, 1] = g[1, 0] = g[2, 3] = g[3, 2] = -0.6
    for i, j in [(0, 2), (2, 1), (1, 3), (3, 0)]:
        g[i, j] = 1.0

    def det_at(x: float) -> complex:
        h = g.copy()
        h[3, 0] = x
        return complex(np.linalg.det(h))

    d0, d1 = det_at(0.0), det_at(1.0)
    g[3, 0] = -d0 / (d1 - d0)
    return LagNetworkToy(g, 1.0)


def ring_plus_pair() -> LagNetworkToy:
    """Ring (w = 0.6) plus an unstable reciprocal pair 0-2: kappa = 2 < d_alg = d_conn = 4."""

    g = ring(w=0.6).g.copy()
    g[0, 2] = g[2, 0] = -1.3
    return LagNetworkToy(g, 1.0)


def scaled_to_boundary(
    toy: LagNetworkToy, target, lo: float = 0.0, hi: float = 10.0
) -> LagNetworkToy:
    """Scale G by c so that portfolio ``target`` sits exactly on its stability boundary."""

    def alpha(c: float) -> float:
        return LagNetworkToy(c * toy.g, toy.a, toy.block).alpha(target)

    if alpha(hi) <= 0:
        raise ValueError("target never becomes unstable in the bracket")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if alpha(mid) > 0:
            hi = mid
        else:
            lo = mid
    return LagNetworkToy(0.5 * (lo + hi) * toy.g, toy.a, toy.block)
