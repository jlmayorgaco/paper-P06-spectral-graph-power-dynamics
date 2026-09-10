"""Spectral composability order: a branch-independent replacement-order measure.

The quantity the earlier work tracked, the real part of one selected mode, needs
a modal label. Labels require matching eigenvectors between systems of different
dimension, which is exactly what broke when the inter-area family acquired two
descendants of indistinguishable shape.

This module defines the replacement order through a COUNT instead:

    N(S)     the number of eigenvalues of the replaced system inside a protected
             spectral region Gamma, with algebraic multiplicity
    kappa    the smallest replacement set whose count differs from the baseline

Nothing in it selects, matches, orders or names an eigenvalue, so nothing in it
can chatter. See ``theory/spectral_composability_order.md`` for the definition
and ``theory/F2B_composability_region_theorem.md`` for what it obeys.

``kappa`` says how many replacements are too many; it does not say which. The
minimal spectral incompatibility hypergraph ``H`` (the minimal sets that change
the count) keeps the coalitions themselves; ``kappa`` is its smallest hyperedge.
See ``theory/F2C_incompatibility_hypergraph.md``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import combinations
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

#: Returned when no tested subset changes the count. Kept as a sentinel rather
#: than ``math.inf`` so that it survives a round trip through CSV as an integer.
NO_FINITE_ORDER = -1


class SpectralRegion(Protocol):
    """A protected region of the complex plane."""

    def contains(self, values: NDArray[np.complex128]) -> NDArray[np.bool_]: ...


@dataclass(frozen=True)
class BandRightHalfPlane:
    """``Gamma``: the open right half plane intersected with a frequency band.

    ``half_plane_only`` keeps one representative per conjugate pair, which is the
    convention used throughout this project. It changes every count by a factor
    of two and changes no verdict, but it must be fixed once and stated.
    """

    low_hz: float
    high_hz: float
    zero_magnitude: float = 1e-3
    half_plane_only: bool = True

    def contains(self, values):
        v = np.asarray(values, dtype=np.complex128)
        frequency = np.abs(v.imag) / (2.0 * np.pi)
        inside = (
            (v.real > 0.0)
            & (frequency >= self.low_hz)
            & (frequency <= self.high_hz)
            & (np.abs(v) > self.zero_magnitude)
        )
        if self.half_plane_only:
            inside &= v.imag >= 0.0
        return inside

    def distance_to_boundary(self, values) -> float:
        """Largest real part among band modes: the signed distance to dGamma.

        Positive means at least one band mode is already inside Gamma. A change of
        sign of this quantity is exactly a crossing of the imaginary-axis part of
        the boundary, which is the transition the region theorem describes.
        """

        v = np.asarray(values, dtype=np.complex128)
        frequency = np.abs(v.imag) / (2.0 * np.pi)
        keep = (
            (frequency >= self.low_hz)
            & (frequency <= self.high_hz)
            & (np.abs(v) > self.zero_magnitude)
        )
        if self.half_plane_only:
            keep &= v.imag >= 0.0
        return float(v.real[keep].max()) if keep.any() else float("-inf")


def spectral_count(values, region: SpectralRegion) -> int:
    """``N_Gamma``: eigenvalues inside the region, with algebraic multiplicity."""

    return int(np.count_nonzero(region.contains(values)))


def composability_order(
    counts: Mapping[frozenset, int], *, baseline: frozenset | None = None
) -> int:
    """``kappa_Gamma``: the smallest set whose count differs from the baseline.

    ``counts`` maps each tested replacement set to its spectral count. Returns
    ``NO_FINITE_ORDER`` when no tested set differs, which is a statement about
    the tested family and never about untested sets.
    """

    empty = frozenset() if baseline is None else baseline
    if empty not in counts:
        raise ValueError("the baseline set must be present in the counts")
    reference = counts[empty]
    sizes = sorted({len(s) for s in counts if s != empty})
    for size in sizes:
        if any(
            counts[s] != reference for s in counts if len(s) == size and s != empty
        ):
            return size
    return NO_FINITE_ORDER


def witnesses(
    counts: Mapping[frozenset, int], *, baseline: frozenset | None = None
) -> tuple[frozenset, ...]:
    """Every minimum-size set that changes the count."""

    order = composability_order(counts, baseline=baseline)
    if order == NO_FINITE_ORDER:
        return ()
    empty = frozenset() if baseline is None else baseline
    reference = counts[empty]
    return tuple(
        sorted(
            (s for s in counts if len(s) == order and counts[s] != reference),
            key=lambda s: sorted(s),
        )
    )


def _check_downward_closed(tested: Iterable[frozenset]) -> None:
    family = set(tested)
    for s in family:
        for member in s:
            if s - {member} not in family:
                raise ValueError(
                    f"tested family is not downward closed: {sorted(s)} is tested "
                    f"but {sorted(s - {member})} is not"
                )


def unsafe_family(
    counts: Mapping[frozenset, int], *, baseline: frozenset | None = None
) -> frozenset[frozenset]:
    """``U``: every tested set whose count differs from the baseline.

    For a baseline-stable region (baseline count zero) this is exactly the family
    of sets with ``N_Gamma(S) > 0``.
    """

    empty = frozenset() if baseline is None else baseline
    if empty not in counts:
        raise ValueError("the baseline set must be present in the counts")
    reference = counts[empty]
    return frozenset(s for s in counts if s != empty and counts[s] != reference)


def incompatibility_hypergraph(
    counts: Mapping[frozenset, int], *, baseline: frozenset | None = None
) -> tuple[frozenset, ...]:
    """``H_Gamma``: the minimal members of the unsafe family.

    A hyperedge is a tested set that disturbs Gamma while none of its proper
    subsets does. The tested family must be downward closed, otherwise "every
    proper subset" would silently range over untested sets. Hyperedges are
    returned in canonical order: by size, then lexicographically.

    ``H`` is an antichain by construction and ``min |S|`` over it is
    ``composability_order``. See ``theory/F2C_incompatibility_hypergraph.md``.
    """

    _check_downward_closed(counts)
    unsafe = unsafe_family(counts, baseline=baseline)
    minimal = [s for s in unsafe if not any(r < s for r in unsafe)]
    return tuple(sorted(minimal, key=lambda s: (len(s), sorted(s))))


def hypergraph_order(edges: Iterable[frozenset]) -> int:
    """``kappa_Gamma`` read off the hypergraph: the smallest hyperedge."""

    sizes = [len(s) for s in edges]
    return min(sizes) if sizes else NO_FINITE_ORDER


def hypergraph_label(edges: Iterable[frozenset]) -> str:
    """Canonical text label, e.g. ``30+33+35|30+37``; ``EMPTY`` when ``H`` is empty."""

    ordered = sorted(edges, key=lambda s: (len(s), sorted(s)))
    if not ordered:
        return "EMPTY"
    return "|".join("+".join(str(m) for m in sorted(s)) for s in ordered)


def parse_hypergraph_label(label: str) -> tuple[frozenset, ...]:
    if label == "EMPTY":
        return ()
    return tuple(
        frozenset(int(m) for m in part.split("+")) for part in label.split("|")
    )


def is_antichain(edges: Iterable[frozenset]) -> bool:
    family = list(edges)
    return not any(a < b for a in family for b in family)


def upward_closure_defect(
    counts: Mapping[frozenset, int], *, baseline: frozenset | None = None
) -> tuple[frozenset, ...]:
    """Tested sets that CONTAIN a hyperedge and are nevertheless safe.

    Empty exactly when the unsafe family is upward closed within the tested
    family, which is when the safe family is a simplicial complex whose minimal
    non-faces are ``H``. Nothing in the dynamics guarantees that, so it is
    measured rather than assumed.
    """

    edges = incompatibility_hypergraph(counts, baseline=baseline)
    unsafe = unsafe_family(counts, baseline=baseline)
    empty = frozenset() if baseline is None else baseline
    return tuple(
        sorted(
            (
                s
                for s in counts
                if s != empty and s not in unsafe and any(e <= s for e in edges)
            ),
            key=lambda s: (len(s), sorted(s)),
        )
    )


def maximal_free_sets(
    edges: Iterable[frozenset], ground: Iterable
) -> tuple[frozenset, ...]:
    """Maximal portfolios containing no hyperedge: guaranteed safe.

    Every unsafe set contains a hyperedge, so a hyperedge-free portfolio is safe
    whether or not stability is monotone. These are the complements of the
    minimal transversals of ``H``: the smallest sets of machines to retain.
    """

    family = [frozenset(e) for e in edges]
    ground = sorted(set(ground))
    free = [
        frozenset(c)
        for k in range(len(ground), -1, -1)
        for c in combinations(ground, k)
        if not any(e <= frozenset(c) for e in family)
    ]
    maximal = [s for s in free if not any(s < t for t in free)]
    return tuple(sorted(maximal, key=lambda s: (-len(s), sorted(s))))


def hypergraph_move(
    before: Iterable[frozenset], after: Iterable[frozenset]
) -> dict[str, object]:
    """Describe how ``H`` changed across one transition.

    ``entered`` and ``left`` are the hyperedges gained and lost. The move is a
    ``CONTRACTION`` when every lost hyperedge strictly contains an entering one
    and nothing else changed, an ``EXPANSION`` for the reverse, and otherwise it
    is described by what entered and left.
    """

    old, new = set(before), set(after)
    entered = sorted(new - old, key=lambda s: (len(s), sorted(s)))
    left = sorted(old - new, key=lambda s: (len(s), sorted(s)))
    if not entered and not left:
        kind = "NONE"
    elif entered and left and all(any(e < s for e in entered) for s in left) and all(
        any(e < s for s in left) for e in entered
    ):
        kind = "CONTRACTION"
    elif entered and left and all(any(s < e for e in entered) for s in left) and all(
        any(s < e for s in left) for e in entered
    ):
        kind = "EXPANSION"
    elif entered and not left:
        kind = "INSERTION"
    elif left and not entered:
        kind = "DELETION"
    else:
        kind = "REORGANISATION"
    return {"kind": kind, "entered": tuple(entered), "left": tuple(left)}


def mobius_terms(values: Mapping[frozenset, float]) -> dict[frozenset, float]:
    """Moebius transform of a set function over the lattice of tested sets.

    Provided so that the composability order and the interaction order can be
    computed from the same data and compared. They are logically independent;
    ``tests/test_spectral_composability.py`` carries counterexamples both ways.
    """

    out: dict[frozenset, float] = {}
    for subset in values:
        total = 0.0
        for other in values:
            if other <= subset:
                total += (-1.0) ** (len(subset) - len(other)) * values[other]
        out[subset] = total
    return out


def interaction_order(values: Mapping[frozenset, float], *, tolerance: float = 0.0) -> int:
    """Largest set size carrying a non-vanishing Moebius term.

    This is the ``irreducible connected interaction order``. It answers "how much
    of the change is genuinely joint", which is NOT the question
    ``composability_order`` answers.
    """

    terms = mobius_terms(values)
    orders = [len(s) for s, v in terms.items() if abs(v) > tolerance and s]
    return max(orders) if orders else 0
