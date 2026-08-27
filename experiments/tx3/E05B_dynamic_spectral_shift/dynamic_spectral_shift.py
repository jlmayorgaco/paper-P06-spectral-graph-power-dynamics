from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.linalg import eigvals, solve
from scipy.sparse import csc_matrix


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
sys.path[:0] = [str(ROOT), str(E03)]

from campaign_core import OperatorPoint, coalition_externality  # noqa: E402
from dicgrid.connected.sparse_logdet import sparse_logabsdet  # noqa: E402


ACTIONS = tuple(range(1, 9))


def all_vertices() -> list[tuple[int, ...]]:
    return [()] + [subset for order in (1, 2, 3) for subset in combinations(ACTIONS, order)]


def all_coalitions() -> list[tuple[int, ...]]:
    return [coalition for order in (2, 3) for coalition in combinations(ACTIONS, order)]


def mobius_vertices(coalition: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [
        tuple(coalition[index] for index in range(len(coalition)) if mask & (1 << index))
        for mask in range(1 << len(coalition))
    ]


def subset_key(subset: Iterable[int]) -> str:
    labels = tuple(subset)
    return "EMPTY" if not labels else "-".join(f"A{value}" for value in labels)


def parse_subset_key(value: str) -> tuple[int, ...]:
    if value in {"", "EMPTY"}:
        return ()
    return tuple(int(label.removeprefix("A")) for label in value.replace(",", "-").split("-"))


def reduced_spectrum(point: OperatorPoint) -> tuple[np.ndarray, float, float]:
    """Return finite poles, log|det(-gy)|, and log|det(M)|.

    This evaluates each factor independently from the full descriptor logdet.
    The finite poles are eigenvalues of M^-1(fx-fy gy^-1 gx).
    """

    if point.status != "SUCCESS" or point.jacobian is None or point.descriptor is None:
        raise ValueError(f"operator point is unavailable: {point.status}")
    n = int(point.metadata["dynamic_state_count"])
    m = int(point.metadata["algebraic_variable_count"])
    jacobian = point.jacobian.toarray()
    fx = jacobian[:n, :n]
    fy = jacobian[:n, n : n + m]
    gx = jacobian[n : n + m, :n]
    gy = jacobian[n : n + m, n : n + m]
    mass = np.asarray(point.descriptor.diagonal()[:n], dtype=float)
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0.0):
        raise FloatingPointError("dynamic descriptor mass is not strictly positive")
    algebraic = sparse_logabsdet(csc_matrix(-gy)).value
    mass_factor = float(np.sum(np.log(mass)))
    algebraic_response = solve(gy, gx, assume_a="gen", check_finite=True)
    reduced = (fx - fy @ algebraic_response) / mass[:, None]
    poles = eigvals(reduced, overwrite_a=True, check_finite=False)
    if not np.isfinite(poles).all():
        raise FloatingPointError("finite-pole spectrum is not finite")
    return np.asarray(poles, dtype=np.complex128), float(algebraic), mass_factor


def pole_logdet(poles: np.ndarray, frequencies_hz: np.ndarray, sigma: float = 0.0) -> np.ndarray:
    s = sigma + 2j * math.pi * np.asarray(frequencies_hz, dtype=float)
    distances = np.abs(s[:, None] - np.asarray(poles, dtype=np.complex128)[None, :])
    if np.any(distances <= 0.0) or not np.isfinite(distances).all():
        raise FloatingPointError("frequency contour intersects a finite pole")
    return np.sum(np.log(distances), axis=1)


def integrated_pole_factor(
    poles: np.ndarray, frequencies_hz: np.ndarray, weights: np.ndarray, sigma: float = 0.0
) -> float:
    return float(np.dot(np.asarray(weights, dtype=float), pole_logdet(poles, frequencies_hz, sigma)))


def component_externality(
    vertex_records: dict[tuple[int, ...], dict[str, float]],
    coalition: tuple[int, ...],
    field: str,
) -> float:
    values = {subset: float(vertex_records[subset][field]) for subset in mobius_vertices(coalition)}
    return coalition_externality(values, coalition)


def connected_resolvent_trace(
    spectra: dict[tuple[int, ...], np.ndarray], coalition: tuple[int, ...], s: complex
) -> complex:
    total = 0.0j
    for subset in mobius_vertices(coalition):
        sign = (-1.0) ** (len(coalition) - len(subset))
        total += sign * np.sum(1.0 / (s - spectra[subset]))
    return complex(total)


def contour_for_reference(baseline: np.ndarray, reference: complex, gap_fraction: float) -> tuple[complex, float]:
    spectrum = np.asarray(baseline, dtype=np.complex128)
    index = int(np.argmin(np.abs(spectrum - reference)))
    center = complex(spectrum[index])
    other = np.delete(spectrum, index)
    gap = float(np.min(np.abs(other - center)))
    if not math.isfinite(gap) or gap <= 0.0:
        raise FloatingPointError("reference pole is not isolated")
    return center, gap_fraction * gap


def exact_contour_moments(
    spectra: dict[tuple[int, ...], np.ndarray], coalition: tuple[int, ...], center: complex, radius: float
) -> dict[str, Any]:
    counts: dict[tuple[int, ...], int] = {}
    sums: dict[tuple[int, ...], complex] = {}
    for subset in mobius_vertices(coalition):
        inside = np.abs(spectra[subset] - center) < radius
        counts[subset] = int(np.count_nonzero(inside))
        sums[subset] = complex(np.sum(spectra[subset][inside]))
    signs = {subset: (-1.0) ** (len(coalition) - len(subset)) for subset in counts}
    mu0 = sum(signs[subset] * counts[subset] for subset in counts)
    mu1 = sum(signs[subset] * sums[subset] for subset in sums)
    return {
        "mu0": complex(mu0),
        "mu1": complex(mu1),
        "vertex_counts": {subset_key(key): value for key, value in counts.items()},
        "single_pole_each_vertex": all(value == 1 for value in counts.values()),
    }


def numerical_contour_moments(
    spectra: dict[tuple[int, ...], np.ndarray],
    coalition: tuple[int, ...],
    center: complex,
    radius: float,
    nodes: int,
) -> tuple[complex, complex]:
    theta = 2.0 * math.pi * np.arange(nodes, dtype=float) / nodes
    unit = np.exp(1j * theta)
    contour = center + radius * unit
    xi = np.asarray([connected_resolvent_trace(spectra, coalition, value) for value in contour])
    mu0 = np.mean(xi * radius * unit)
    mu1 = np.mean(contour * xi * radius * unit)
    return complex(mu0), complex(mu1)


def positive_oscillatory_poles(poles: np.ndarray, band_hz: tuple[float, float]) -> np.ndarray:
    values = np.asarray(poles, dtype=np.complex128)
    frequencies = values.imag / (2.0 * math.pi)
    selected = values[(frequencies >= band_hz[0]) & (frequencies <= band_hz[1])]
    return selected[np.argsort(selected.imag)]
