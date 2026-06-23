"""IEEE 39 validation harness for the rational graph-filter damping estimator.

This script is intentionally conservative. It does not fabricate an IEEE 39
case and it does not claim validation if ANDES or the dynamic case is missing.

Usage after installing ANDES and adding case files:

    python run_ieee39_validation.py --raw cases/ieee39.raw --dyr cases/ieee39.dyr

The functions for the spectral estimator are implemented for a supplied
second-order surrogate (M, D, L). The ANDES extraction layer is left explicit
because the mapping from full IBR DAEs to scalar M,D,L must be documented for
each model.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass
class EstimatorResult:
    zeta_diag: float
    zeta_second_order: float
    corrected_pole_real: float
    corrected_pole_imag: float
    critical_mode: int
    trigger_eta: float
    trigger_reject: bool


def damping_ratio(s: complex) -> float:
    """Return damping ratio -Re(s)/|s| for one pole."""
    mag = abs(s)
    if mag == 0:
        return np.inf
    return float(-s.real / mag)


def oscillatory_margin(poles: Iterable[complex], imag_tol: float = 1e-8) -> float:
    """Minimum damping ratio over oscillatory poles."""
    ratios = [damping_ratio(s) for s in poles if abs(s.imag) > imag_tol and abs(s) > imag_tol]
    return min(ratios) if ratios else np.inf


def second_order_surrogate_estimator(
    M: np.ndarray,
    D: np.ndarray,
    L: np.ndarray,
    eta_threshold: float = 0.25,
    eps: float = 1e-9,
) -> EstimatorResult:
    """Compute diagonal and second-order damping-margin estimates.

    Parameters
    ----------
    M, D, L:
        Second-order surrogate matrices in
        M theta_ddot + D theta_dot + L theta = 0.
    eta_threshold:
        Conservative coupling/degeneracy trigger threshold.

    Notes
    -----
    This estimator assumes symmetric L, diagonal positive M, and a physically
    meaningful local damping matrix D. It is a screening model, not the full
    ANDES DAE.
    """
    M = np.asarray(M, dtype=float)
    D = np.asarray(D, dtype=float)
    L = np.asarray(L, dtype=float)

    if M.ndim == 1:
        Minv_half = np.diag(1.0 / np.sqrt(M))
    else:
        Minv_half = np.diag(1.0 / np.sqrt(np.diag(M)))

    Lt = Minv_half @ L @ Minv_half
    Dt = Minv_half @ D @ Minv_half

    nu, Q = np.linalg.eigh(Lt)
    # Drop the near-zero angle-reference mode.
    valid = np.where(nu > 1e-10)[0]
    if len(valid) == 0:
        raise ValueError("No nonzero inertial-Laplacian modes found.")

    Gamma = Q.T @ Dt @ Q
    delta = np.diag(Gamma)
    E = Gamma - np.diag(delta)

    zeta_diag_modes = np.full_like(nu, np.inf, dtype=float)
    zeta_diag_modes[valid] = delta[valid] / (2.0 * np.sqrt(nu[valid]))
    c = int(np.argmin(zeta_diag_modes))
    zeta_diag = float(zeta_diag_modes[c])

    disc = delta[c] ** 2 - 4.0 * nu[c]
    if disc >= 0:
        # Non-oscillatory diagonal critical candidate; reject cheap estimate.
        s0 = complex((-delta[c] - np.sqrt(disc)) / 2.0, 0.0)
        trigger_reject = True
    else:
        s0 = complex(-delta[c] / 2.0, np.sqrt(4.0 * nu[c] - delta[c] ** 2) / 2.0)
        trigger_reject = False

    correction_sum = 0.0 + 0.0j
    eta = 0.0
    for ell in valid:
        ell = int(ell)
        if ell == c:
            continue
        p_ell = s0 * s0 + delta[ell] * s0 + nu[ell]
        if abs(p_ell) < eps:
            trigger_reject = True
            continue
        correction_sum += E[c, ell] * E[ell, c] / p_ell

        gap = abs(np.sqrt(nu[ell]) - np.sqrt(nu[c]))
        eta = max(eta, abs(E[c, ell]) / (gap + eps))

    ds = (s0 * s0 / (2.0 * s0 + delta[c])) * correction_sum
    s2 = s0 + ds
    zeta_second = damping_ratio(s2)
    trigger_reject = bool(trigger_reject or eta > eta_threshold or abs(ds) > 0.25 * max(abs(s0), eps))

    return EstimatorResult(
        zeta_diag=zeta_diag,
        zeta_second_order=float(zeta_second),
        corrected_pole_real=float(s2.real),
        corrected_pole_imag=float(s2.imag),
        critical_mode=c,
        trigger_eta=float(eta),
        trigger_reject=trigger_reject,
    )


def load_andes_case(raw: Path, dyr: Path):
    """Load an ANDES case if ANDES is installed.

    This function deliberately returns the ANDES system object only. Extraction
    of M,D,L and full poles must be implemented with documented model-specific
    choices.
    """
    try:
        import andes  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "ANDES is not installed. Install it in the active Python environment "
            "before running IEEE 39 validation."
        ) from exc

    if not raw.exists():
        raise FileNotFoundError(raw)
    if not dyr.exists():
        raise FileNotFoundError(dyr)

    ss = andes.load(str(raw), addfile=str(dyr), setup=False)
    ss.setup()
    ss.PFlow.run()
    return ss


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=False, help="IEEE 39 RAW file")
    parser.add_argument("--dyr", type=Path, required=False, help="IEEE 39 DYR file")
    parser.add_argument("--out", type=Path, default=Path("outputs") / "ieee39_pending")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    status = {
        "status": "not_run",
        "reason": "This harness requires ANDES and an IEEE 39 dynamic case.",
        "raw": str(args.raw) if args.raw else None,
        "dyr": str(args.dyr) if args.dyr else None,
        "required_next_step": "Implement model-specific extraction of full poles and M,D,L surrogate.",
    }

    if args.raw and args.dyr:
        try:
            _ = load_andes_case(args.raw, args.dyr)
            status["status"] = "case_loaded"
            status["reason"] = "ANDES case loaded; extraction layer still required."
        except Exception as exc:  # noqa: BLE001 - record validation blocker.
            status["status"] = "blocked"
            status["reason"] = str(exc)

    with (args.out / "validation_status.json").open("w", encoding="utf-8") as f:
        json.dump(status, f, indent=2)

    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
