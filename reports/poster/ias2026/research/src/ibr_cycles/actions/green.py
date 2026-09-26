from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .low_rank import ActionFactor

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class ActionGreen:
    """The Action Green operator K(s) = V^H T_0(s)^-1 U of a portfolio.

    Block K_ab(s) reads: action b injects into the controller-dressed grid and
    is observed by action a. It is the only place where the network and the
    hidden controller dynamics enter the interaction calculus.
    """

    a0: ComplexMatrix
    factors: tuple[ActionFactor, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(factor.name for factor in self.factors)

    @property
    def block_slices(self) -> tuple[slice, ...]:
        offsets = np.cumsum([0] + [factor.rank for factor in self.factors])
        return tuple(
            slice(int(offsets[i]), int(offsets[i + 1]))
            for i in range(len(self.factors))
        )

    @property
    def dimension(self) -> int:
        return int(sum(factor.rank for factor in self.factors))

    def stacked(self) -> tuple[ComplexMatrix, ComplexMatrix, ComplexMatrix]:
        u = np.hstack([factor.u for factor in self.factors])
        v = np.hstack([factor.v for factor in self.factors])
        c = np.zeros((self.dimension, self.dimension), dtype=np.complex128)
        for factor, block in zip(self.factors, self.block_slices, strict=True):
            c[block, block] = factor.c
        return u, c, v

    def t0(self, s: complex) -> ComplexMatrix:
        n = self.a0.shape[0]
        return s * np.eye(n, dtype=np.complex128) - self.a0

    def k(self, s: complex) -> ComplexMatrix:
        u, _, v = self.stacked()
        return v.conj().T @ np.linalg.solve(self.t0(s), u)

    def m(self, s: complex) -> ComplexMatrix:
        """M(s) = C K(s), so that det(T_S)/det(T_0) equals det(I + M(s))."""

        _, c, _ = self.stacked()
        return c @ self.k(s)

    def block(self, matrix: ComplexMatrix, a: int, b: int) -> ComplexMatrix:
        slices = self.block_slices
        return matrix[slices[a], slices[b]]

    def determinant_ratio(self, s: complex) -> complex:
        identity = np.eye(self.dimension, dtype=np.complex128)
        return complex(np.linalg.det(identity + self.m(s)))

    def exact_determinant_ratio(self, s: complex) -> complex:
        """Reference value computed without the determinant lemma."""

        delta = np.zeros_like(self.a0)
        for factor in self.factors:
            delta = delta + factor.delta_t()
        numerator = np.linalg.det(self.t0(s) + delta)
        denominator = np.linalg.det(self.t0(s))
        return complex(numerator / denominator)

    def lemma_residual(self, points: Sequence[complex]) -> float:
        """EXACT IDENTITY check of the matrix determinant lemma."""

        return max(
            abs(self.determinant_ratio(s) - self.exact_determinant_ratio(s))
            for s in points
        )

    def subset(self, names: Sequence[str]) -> ActionGreen:
        wanted = list(names)
        chosen = tuple(f for f in self.factors if f.name in wanted)
        if len(chosen) != len(wanted):
            missing = set(wanted) - {f.name for f in chosen}
            raise KeyError(f"unknown actions: {sorted(missing)}")
        return ActionGreen(a0=self.a0, factors=chosen)
