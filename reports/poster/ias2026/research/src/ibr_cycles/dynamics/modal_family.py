"""Inter-area modal FAMILY tracking, and the envelope observable built on it.

Why this module exists. A single-branch tracker picks, out of the post-change
spectrum, the one mode with the highest modal assurance to a frozen reference
shape. That is well posed only while the reference has exactly one descendant.
On the IEEE-39 replacement it does not: the base inter-area mode acquires two
descendants of nearly equal assurance, the selector chatters between them, and
the resulting scalar is discontinuous even though the physics is not
(``docs/V2B_RESULTS.md``, ``FAILED_EXPERIMENTS.md`` F7).

The fix is to stop selecting. The tracked object is the whole descendant set

    C_IA(theta) = { modes in the frozen band whose shape overlaps the frozen
                    base inter-area reference by at least ``threshold`` }

and the stability observable is its envelope

    alpha_IA(theta) = max over lambda in C_IA(theta) of Re(lambda).

This reduces to the old observable whenever the family has one member, so it is
a strict generalization rather than a different quantity.

Identifiability is checked at the level of the SUBSPACE the family spans, not at
the level of its individual members, because individual assurance is exactly what
fails here. ``subspace_similarity`` reports the smallest principal-angle cosine
between two families, which stays near one across a swap of the members.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .modal_tracking import modal_assurance
from .modes import Mode, Spectrum

ComplexMatrix = NDArray[np.complex128]

ZERO_MODE = 1e-3
"""Below this magnitude a mode is one of the reference double zeros, not dynamics."""


def orthonormal(vectors: ComplexMatrix) -> ComplexMatrix:
    """Orthonormal basis of the column span, robust to near-parallel columns."""

    if vectors.size == 0:
        return vectors.reshape(vectors.shape[0], 0)
    u, singular, _ = np.linalg.svd(vectors, full_matrices=False)
    keep = singular > max(vectors.shape) * np.finfo(float).eps * singular[0]
    return u[:, keep]


def principal_cosines(a: ComplexMatrix, b: ComplexMatrix) -> NDArray[np.float64]:
    """Cosines of the principal angles between two column spans, descending."""

    qa, qb = orthonormal(a), orthonormal(b)
    if qa.shape[1] == 0 or qb.shape[1] == 0:
        return np.zeros(0)
    singular = np.linalg.svd(qa.conj().T @ qb, compute_uv=False)
    return np.clip(np.asarray(singular, dtype=float), 0.0, 1.0)


def subspace_similarity(a: ComplexMatrix, b: ComplexMatrix) -> float:
    """Worst principal-angle cosine, over ``min(dim a, dim b)`` directions.

    One means the smaller span sits inside the larger. Unlike a modal assurance
    between individual vectors, this is invariant to any change of basis inside
    either family, so a swap between two descendants leaves it unchanged.
    """

    cosines = principal_cosines(a, b)
    return float(cosines.min()) if cosines.size else 0.0


@dataclass(frozen=True)
class ModalFamily:
    """The descendant set of the frozen inter-area reference, and its envelope."""

    modes: tuple[Mode, ...]
    overlaps: tuple[float, ...]
    basis: ComplexMatrix
    independence: float
    """Smallest singular value of the normalized member vectors stacked as columns.

    Near zero means the members are nearly parallel and the family, while well
    defined as a set, spans less than its cardinality suggests.
    """

    @property
    def size(self) -> int:
        return len(self.modes)

    @property
    def worst(self) -> Mode:
        return max(self.modes, key=lambda m: m.real)

    @property
    def alpha(self) -> float:
        """``alpha_IA``: the envelope of the family, its least damped member."""

        return float(self.worst.real)

    @property
    def frequency_worst_hz(self) -> float:
        return float(self.worst.frequency_hz)

    @property
    def best(self) -> Mode:
        return min(self.modes, key=lambda m: m.real)

    @property
    def frequencies_hz(self) -> tuple[float, ...]:
        return tuple(m.frequency_hz for m in sorted(self.modes, key=lambda x: x.real))

    @property
    def reals(self) -> tuple[float, ...]:
        return tuple(m.real for m in sorted(self.modes, key=lambda x: x.real))

    @property
    def max_condition(self) -> float:
        return float(max(m.condition for m in self.modes))


def band_candidates(
    spectrum: Spectrum, band_hz: tuple[float, float]
) -> tuple[Mode, ...]:
    """Modes inside the frozen band, one representative per conjugate pair."""

    return tuple(
        mode
        for mode in spectrum.modes
        if abs(mode.value) > ZERO_MODE
        and mode.value.imag >= 0.0
        and band_hz[0] <= mode.frequency_hz <= band_hz[1]
    )


def restrict(mode: Mode, positions: list[int]) -> NDArray[np.complex128]:
    vector = mode.right[positions]
    norm = float(np.linalg.norm(vector))
    return vector / norm if norm > 0.0 else vector


def descendant_family(
    reference: NDArray[np.complex128],
    positions: list[int],
    spectrum: Spectrum,
    *,
    band_hz: tuple[float, float],
    threshold: float,
) -> ModalFamily | None:
    """Every band mode whose shape overlaps the frozen reference by ``threshold``.

    ``reference`` is the frozen inter-area shape restricted to the coordinates
    shared by the reference and the candidate case; ``positions`` locates those
    same coordinates in the candidate. Returns ``None`` when the family is empty,
    which the caller must classify as a tracking failure rather than fall back to
    the nearest frequency.
    """

    kept: list[Mode] = []
    overlaps: list[float] = []
    for mode in band_candidates(spectrum, band_hz):
        overlap = modal_assurance(reference, mode.right[positions])
        if overlap >= threshold:
            kept.append(mode)
            overlaps.append(float(overlap))
    if not kept:
        return None
    columns = np.column_stack([restrict(mode, positions) for mode in kept])
    singular = np.linalg.svd(columns, compute_uv=False)
    return ModalFamily(
        modes=tuple(kept),
        overlaps=tuple(overlaps),
        basis=orthonormal(columns),
        independence=float(singular[-1]),
    )
