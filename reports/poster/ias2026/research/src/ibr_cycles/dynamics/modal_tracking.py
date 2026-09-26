from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linear_sum_assignment

from .modes import Mode, Spectrum


def modal_assurance(
    left: NDArray[np.complex128], right: NDArray[np.complex128]
) -> float:
    """Modal assurance criterion between two eigenvectors (1.0 means identical)."""

    numerator = abs(complex(left.conj() @ right)) ** 2
    denominator = float((left.conj() @ left).real * (right.conj() @ right).real)
    if denominator <= 0.0:
        return 0.0
    return float(numerator / denominator)


def track(
    reference: Spectrum,
    candidate: Spectrum,
    *,
    distance_weight: float = 1.0,
    mac_weight: float = 1.0,
) -> tuple[tuple[int, ...], tuple[float, ...]]:
    """Match candidate modes to reference modes by distance and eigenvector overlap.

    Sorting by frequency alone is never used because it reorders coalescing
    modes. Returns the assignment and the MAC of each accepted pair, so that a
    caller can reject a weak track instead of reporting a wrong one.
    """

    n = len(reference.modes)
    m = len(candidate.modes)
    cost = np.zeros((n, m))
    scale = max(float(np.abs(reference.values).max()), 1.0)
    for i, ref in enumerate(reference.modes):
        for j, cand in enumerate(candidate.modes):
            gap = abs(ref.value - cand.value) / scale
            mac = modal_assurance(ref.right, cand.right)
            cost[i, j] = distance_weight * gap + mac_weight * (1.0 - mac)
    rows, cols = linear_sum_assignment(cost)
    assignment = [-1] * n
    quality = [0.0] * n
    for i, j in zip(rows, cols, strict=True):
        assignment[i] = int(j)
        quality[i] = modal_assurance(reference.modes[i].right, candidate.modes[j].right)
    return tuple(assignment), tuple(quality)


def coalescence_warning(
    spectrum: Spectrum, *, condition_limit: float = 1e6
) -> list[Mode]:
    """Modes whose eigenvalue condition number makes tracking untrustworthy."""

    return [mode for mode in spectrum.modes if mode.condition > condition_limit]
