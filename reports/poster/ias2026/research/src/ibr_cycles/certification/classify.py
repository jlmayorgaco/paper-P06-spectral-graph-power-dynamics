"""Physical stability classification with frozen numerical criteria (BC00-B).

Four states: STABLE, UNSTABLE, BOUNDARY_OR_UNRESOLVED, INFEASIBLE.

No eigenvalue is removed because it is small. With eps the size of the matrix
error (finite-difference Jacobian error plus the backward error of the
eigensolver) and SAFETY a frozen factor:

Resolved count
    If the distance to the imaginary axis exceeds the error,

        d_axis = min_w sigma_min(A - i w I) > SAFETY * eps,

    then no perturbation of that size moves an eigenvalue across the axis (the
    classical complex stability radius). The number of right-half-plane
    eigenvalues is therefore exact: STABLE if it is 0, UNSTABLE otherwise.
    d_axis is computed numerically: the imaginary part of every eigenvalue, 0,
    every eigenvalue with Re > -0.5, and a 24-point log grid,
    then bounded refinement around the two smallest values.
    It is a numerical criterion, not an interval certificate.
Otherwise
    UNSTABLE if some eigenvalue is resolved positive to first order:
    Re(lambda) > tau with tau = SAFETY * kappa(lambda) * eps and
    kappa = ||x|| ||y|| / |y^H x|. The count is then a lower bound. A defective
    eigenvalue has kappa = inf and is never resolved by this route.
BOUNDARY_OR_UNRESOLVED
    Everything else: an exact zero (Jordan or neutral mode), an eigenvalue
    within the error of the axis, or a non-normal near-axis pseudospectrum.

SAFETY = 10 was fixed in configs/binary_certification_v1/BC00_config.yaml
before any benchmark hypergraph was recomputed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eig, matrix_balance
from scipy.optimize import minimize_scalar

STABLE = "STABLE"
UNSTABLE = "UNSTABLE"
BOUNDARY = "BOUNDARY_OR_UNRESOLVED"
INFEASIBLE = "INFEASIBLE"
SAFETY = 10.0
_EPS = float(np.finfo(np.float64).eps)
WINDOWS = (1e-2, 1e-3, 1e-4, 1e-6, 1e-8)


@dataclass(frozen=True)
class EigenRecord:
    value: complex
    kappa: float
    tau: float
    sign: str  # POSITIVE (resolved), NEGATIVE (Re < 0), NEAR_AXIS


@dataclass(frozen=True)
class SpectrumReport:
    status: str
    n_positive: int  # resolved positive, with multiplicity (a pair counts 2)
    n_positive_real: int
    n_rhp_raw: int  # plain count Re > 0, no tolerance (diagnostic only)
    abscissa: float
    d_axis: float
    threshold: float
    records: tuple[EigenRecord, ...] = field(repr=False)

    def near_axis(self) -> list[EigenRecord]:
        return [r for r in self.records if r.sign == "NEAR_AXIS"]


def _sigma_min(a: NDArray, w: float) -> float:
    n = a.shape[0]
    return float(np.linalg.svd(a - 1j * w * np.eye(n), compute_uv=False)[-1])


def distance_to_axis(a: NDArray, values: NDArray) -> float:
    """min_w sigma_min(A - i w I), numerically (see module docstring)."""

    top = float(np.abs(values.imag).max()) if values.size else 1.0
    slow = values[values.real > -0.5]
    grid = np.unique(
        np.concatenate(
            [
                [0.0],
                np.abs(slow.imag),
                np.geomspace(1e-3, max(1.5 * top, 1.0), 24),
            ]
        )
    )
    sig = np.array([_sigma_min(a, w) for w in grid])
    best = float(sig.min())
    for k in np.argsort(sig)[:2]:
        lo = grid[max(k - 1, 0)]
        hi = grid[min(k + 1, grid.size - 1)]
        if hi > lo:
            res = minimize_scalar(
                lambda w: _sigma_min(a, w),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-6 * max(1.0, hi)},
            )
            best = min(best, float(res.fun))
    return best


def classify_spectrum(
    a: NDArray[np.float64],
    matrix_error: float | NDArray[np.float64] = 0.0,
    safety: float = SAFETY,
    known_zero: bool = False,
) -> SpectrumReport:
    """Classify A. ``matrix_error`` is a norm or, better, the error MATRIX dA.

    A is first balanced by a diagonal similarity D (powers of two). The error
    matrix is transformed with it, D^-1 dA D, so that the error and the
    distance to the axis are measured in the same metric. A scalar error cannot
    be transformed and is used as given.

    ``known_zero``: the caller has an exact zero eigenvector (the verified
    neutral-frequency mode). d_axis is then 0 by construction and is not
    computed; the status is BOUNDARY_OR_UNRESOLVED unless an eigenvalue is
    resolved positive.
    """

    a = np.asarray(a, dtype=np.float64)
    if a.size == 0:
        return SpectrumReport(STABLE, 0, 0, 0, -np.inf, np.inf, 0.0, ())
    a, (scale, _perm) = matrix_balance(a, permute=False, separate=True)
    if np.ndim(matrix_error) == 2:
        err = np.asarray(matrix_error) / scale[:, None] * scale[None, :]
        matrix_error = float(np.linalg.norm(err, 2))
    values, left, right = eig(a, left=True, right=True)
    norm = max(float(np.linalg.norm(a, 2)), 1e-300)
    residuals = np.linalg.norm(a @ right - right * values, axis=0) / (
        norm * np.linalg.norm(right, axis=0)
    )
    eps = float(matrix_error + max(residuals.max(), _EPS) * norm)
    threshold = safety * eps
    records = []
    for k, lam in enumerate(values):
        x, y = right[:, k], left[:, k]
        overlap = abs(np.vdot(y, x))
        kappa = (
            np.inf
            if overlap < 1e-300
            else float(np.linalg.norm(x) * np.linalg.norm(y) / overlap)
        )
        tau = safety * kappa * eps
        if lam.real > tau:
            sign = "POSITIVE"
        elif lam.real < 0 and lam.real < -tau:
            sign = "NEGATIVE"
        else:
            sign = "NEAR_AXIS"
        records.append(EigenRecord(complex(lam), kappa, float(tau), sign))
    d_axis = 0.0 if known_zero else distance_to_axis(a, values)
    if d_axis > threshold:
        # the RHP count is invariant under every admissible perturbation:
        # nothing is near the axis, every sign is resolved
        records = [
            EigenRecord(
                r.value, r.kappa, r.tau, "POSITIVE" if r.value.real > 0 else "NEGATIVE"
            )
            for r in records
        ]
    n_pos = sum(r.sign == "POSITIVE" for r in records)
    n_pos_real = sum(
        r.sign == "POSITIVE" and abs(r.value.imag) <= max(r.tau, 1e-12) for r in records
    )
    if n_pos:
        status = UNSTABLE
    elif d_axis > threshold:
        status = STABLE
    else:
        status = BOUNDARY
    return SpectrumReport(
        status=status,
        n_positive=int(n_pos),
        n_positive_real=int(n_pos_real),
        n_rhp_raw=int(np.count_nonzero(values.real > 0)),
        abscissa=float(values.real.max()),
        d_axis=float(d_axis),
        threshold=float(threshold),
        records=tuple(records),
    )


def near_zero_table(values: NDArray[np.complex128]) -> dict[str, int]:
    """How many eigenvalues fall in each diagnostic window (none is removed)."""

    mags = np.abs(np.asarray(values))
    return {f"|lambda|<={w:g}": int(np.count_nonzero(mags <= w)) for w in WINDOWS}
