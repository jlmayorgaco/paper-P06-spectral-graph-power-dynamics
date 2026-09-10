from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np


def _mobius_sum(values: dict[frozenset[str], float], members: tuple[str, ...]) -> float:
    total = 0.0
    for size in range(len(members) + 1):
        sign = (-1.0) ** (len(members) - size)
        for subset in combinations(members, size):
            total += sign * values[frozenset(subset)]
    return total


@dataclass(frozen=True)
class BranchTrace:
    """Homotopy continuation record of a branch-consistent complex logarithm."""

    amplitudes: tuple[float, ...]
    values: tuple[complex, ...]
    max_step_angle: float
    unwrapped: complex

    @property
    def reliable(self) -> bool:
        """A step that turns by more than pi/2 may have skipped a branch."""

        return self.max_step_angle < np.pi / 2.0


def continued_log(
    ratio: Callable[[float], complex],
    *,
    steps: int = 256,
) -> BranchTrace:
    """Unwrap log of a determinant ratio by continuation from amplitude 0 to 1.

    Independent principal-branch logarithms are wrong whenever the ratio winds
    around the origin, which is precisely the regime of interest, so the branch
    is followed rather than assumed.
    """

    amplitudes = np.linspace(0.0, 1.0, steps + 1)
    values = np.array([ratio(float(a)) for a in amplitudes], dtype=np.complex128)
    if np.any(np.abs(values) < 1e-300):
        raise ValueError("determinant ratio vanishes on the continuation path")
    increments = np.angle(values[1:] / values[:-1])
    phase = float(np.angle(values[0])) + float(np.sum(increments))
    unwrapped = complex(np.log(abs(values[-1])), phase)
    return BranchTrace(
        amplitudes=tuple(float(a) for a in amplitudes),
        values=tuple(complex(v) for v in values),
        max_step_angle=float(np.abs(increments).max()),
        unwrapped=unwrapped,
    )


@dataclass(frozen=True)
class MobiusDecomposition:
    """Finite-amplitude interaction decomposition over the Boolean lattice.

        mu(S) = sum over R contained in S of (-1)^(abs(S)-abs(R)) v(R)

    Order one is the isolated effect, order two the irreducible pair
    interaction, order three the genuine finite-amplitude triple term. A
    nonzero third-order term does not by itself prove that the pure third-order
    term caused instability; compare the truncated reconstructions for that.
    """

    values: dict[frozenset[str], float]
    terms: dict[frozenset[str], float]
    labels: tuple[str, ...]

    def order(self, k: int) -> dict[frozenset[str], float]:
        return {key: value for key, value in self.terms.items() if len(key) == k}

    def truncated(self, members: tuple[str, ...], k: int) -> float:
        """Reconstruction of v(members) keeping interaction orders up to k."""

        total = 0.0
        for size in range(k + 1):
            for subset in combinations(members, size):
                total += self.terms[frozenset(subset)]
        return total

    def reconstruction_error(self, members: tuple[str, ...]) -> float:
        full = self.truncated(members, len(members))
        return float(abs(full - self.values[frozenset(members)]))


def decompose(
    functional: Callable[[tuple[str, ...]], float], labels: Sequence[str]
) -> MobiusDecomposition:
    """Evaluate the functional on every subset and invert the Boolean lattice."""

    names = tuple(labels)
    values: dict[frozenset[str], float] = {}
    for size in range(len(names) + 1):
        for subset in combinations(names, size):
            values[frozenset(subset)] = float(functional(subset))
    terms: dict[frozenset[str], float] = {}
    for size in range(len(names) + 1):
        for subset in combinations(names, size):
            terms[frozenset(subset)] = _mobius_sum(values, subset)
    return MobiusDecomposition(values=values, terms=terms, labels=names)
