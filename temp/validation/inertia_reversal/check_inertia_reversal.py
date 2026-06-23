"""Check whether adding inertia can reduce damping margin.

The script separates two mechanisms:

1. fixed_physical_damping: M_i is increased while D_i is held fixed. This can
   reduce damping ratio even for a scalar oscillator because D_i/M_i falls; it
   is not a non-proportional modal-rotation result.
2. fixed_damping_per_inertia: M_i is increased while d_i = D_i/M_i is held
   fixed. The direct damping dilution is removed. A negative derivative here is
   a genuine modal-rotation/non-proportional effect in the diagonal modal
   screen; the full QEP margin is checked by finite differences.

Outputs a reproducible example and Monte Carlo counts.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "inertia_reversal"


@dataclass
class ReversalExample:
    trial: int
    bus: int
    n: int
    assumption: str
    base_full_zeta: float
    post_full_zeta: float
    full_fd_derivative: float
    base_diag_zeta: float
    post_diag_zeta: float
    diag_fd_derivative: float
    diag_analytic_derivative: float
    critical_mode: int
    min_gap: float
    m_i: float
    d_over_m_i: float


def connected_laplacian(rng: np.random.Generator, n: int) -> np.ndarray:
    """Random connected weighted Laplacian."""
    W = np.zeros((n, n))
    # random tree
    for j in range(1, n):
        i = int(rng.integers(0, j))
        w = float(rng.uniform(0.5, 4.0))
        W[i, j] = W[j, i] = w
    # extra edges
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] == 0 and rng.random() < 0.35:
                W[i, j] = W[j, i] = float(rng.uniform(0.2, 3.0))
    return np.diag(W.sum(axis=1)) - W


def normalized_matrices(L: np.ndarray, M: np.ndarray, D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    Sinv = np.diag(1.0 / np.sqrt(M))
    Lt = Sinv @ L @ Sinv
    Dt = Sinv @ np.diag(D) @ Sinv
    return Lt, Dt


def modal_data(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    Lt, Dt = normalized_matrices(L, M, D)
    vals, Q = np.linalg.eigh(Lt)
    keep = vals > 1e-9
    vals = vals[keep]
    Q = Q[:, keep]
    Gamma = Q.T @ Dt @ Q
    delta = np.diag(Gamma)
    zeta = delta / (2.0 * np.sqrt(vals))
    return Lt, Dt, vals, Q, Gamma, zeta


def full_qep_margin(L: np.ndarray, M: np.ndarray, D: np.ndarray) -> float:
    Lt, Dt = normalized_matrices(L, M, D)
    n = L.shape[0]
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Lt, -Dt]])
    ev = np.linalg.eigvals(A)
    zetas = [-s.real / abs(s) for s in ev if abs(s.imag) > 1e-7 and abs(s) > 1e-9 and s.real < 1e-8]
    return float(min(zetas)) if zetas else 1.0


def dlt_dmi(Lt: np.ndarray, M_i: float, bus: int) -> np.ndarray:
    E = np.zeros_like(Lt)
    E[bus, bus] = 1.0
    return -(E @ Lt + Lt @ E) / (2.0 * M_i)


def analytic_diag_derivative_fixed_ratio(
    L: np.ndarray, M: np.ndarray, dbar: np.ndarray, bus: int, mode: int
) -> float:
    """d zeta_k / d M_i with dbar=D/M fixed, for a simple modal eigenvalue."""
    D = M * dbar
    Lt, Dt, nu, Q, Gamma, zeta = modal_data(L, M, D)
    k = mode
    dLt = dlt_dmi(Lt, M[bus], bus)
    # dDt is zero because D_i/M_i is held fixed.
    dnu = float(Q[:, k].T @ dLt @ Q[:, k])
    ddelta = 0.0
    for ell in range(len(nu)):
        if ell == k:
            continue
        gap = nu[k] - nu[ell]
        if abs(gap) < 1e-10:
            continue
        rot = float(Q[:, ell].T @ dLt @ Q[:, k])
        damp_couple = float(Q[:, ell].T @ Dt @ Q[:, k])
        ddelta += 2.0 * rot * damp_couple / gap
    return float(ddelta / (2.0 * np.sqrt(nu[k])) - zeta[k] * dnu / (2.0 * nu[k]))


def matched_mode_zeta(L: np.ndarray, M: np.ndarray, D: np.ndarray, q_ref: np.ndarray) -> tuple[float, int, float]:
    _, _, nu, Q, _, zeta = modal_data(L, M, D)
    dots = np.abs(Q.T @ q_ref)
    idx = int(np.argmax(dots))
    return float(zeta[idx]), idx, float(dots[idx])


def evaluate_case(
    L: np.ndarray,
    M: np.ndarray,
    dbar: np.ndarray,
    bus: int,
    assumption: str,
    eps: float = 1e-4,
) -> ReversalExample | None:
    if assumption == "fixed_damping_per_inertia":
        D = M * dbar
        M2 = M.copy()
        M2[bus] *= 1.0 + eps
        D2 = M2 * dbar
    elif assumption == "fixed_physical_damping":
        D = M * dbar
        M2 = M.copy()
        M2[bus] *= 1.0 + eps
        D2 = D.copy()
    else:
        raise ValueError(assumption)

    _, _, nu, Q, _, zeta = modal_data(L, M, D)
    mode = int(np.argmin(zeta))
    z0_diag = float(zeta[mode])
    z1_diag, _, overlap = matched_mode_zeta(L, M2, D2, Q[:, mode])
    if overlap < 0.98:
        return None
    fd_diag = (z1_diag - z0_diag) / (M2[bus] - M[bus])

    if assumption == "fixed_damping_per_inertia":
        analytic = analytic_diag_derivative_fixed_ratio(L, M, dbar, bus, mode)
    else:
        analytic = float("nan")

    z0_full = full_qep_margin(L, M, D)
    z1_full = full_qep_margin(L, M2, D2)
    fd_full = (z1_full - z0_full) / (M2[bus] - M[bus])
    gaps = np.diff(np.sort(nu))
    min_gap = float(np.min(gaps)) if len(gaps) else float("nan")

    return ReversalExample(
        trial=-1,
        bus=bus,
        n=L.shape[0],
        assumption=assumption,
        base_full_zeta=z0_full,
        post_full_zeta=z1_full,
        full_fd_derivative=float(fd_full),
        base_diag_zeta=z0_diag,
        post_diag_zeta=z1_diag,
        diag_fd_derivative=float(fd_diag),
        diag_analytic_derivative=float(analytic),
        critical_mode=mode,
        min_gap=min_gap,
        m_i=float(M[bus]),
        d_over_m_i=float(dbar[bus]),
    )


def main() -> None:
    rng = np.random.default_rng(20260608)
    OUT.mkdir(parents=True, exist_ok=True)

    examples: list[ReversalExample] = []
    counts = {
        "fixed_damping_per_inertia": {"tested": 0, "diag_reversal": 0, "full_reversal": 0, "both": 0},
        "fixed_physical_damping": {"tested": 0, "diag_reversal": 0, "full_reversal": 0, "both": 0},
    }

    for trial in range(2500):
        n = int(rng.integers(4, 10))
        L = connected_laplacian(rng, n)
        M = rng.uniform(2.0, 12.0, n)
        # deliberately heterogeneous local damping per inertia
        dbar = rng.lognormal(mean=-0.15, sigma=0.75, size=n)
        # avoid extremely damped/undamped cases
        dbar = np.clip(dbar, 0.08, 2.0)

        for assumption in counts:
            for bus in range(n):
                ex = evaluate_case(L, M, dbar, bus, assumption)
                if ex is None:
                    continue
                ex.trial = trial
                counts[assumption]["tested"] += 1
                diag_rev = ex.diag_fd_derivative < -1e-7
                full_rev = ex.full_fd_derivative < -1e-7
                counts[assumption]["diag_reversal"] += int(diag_rev)
                counts[assumption]["full_reversal"] += int(full_rev)
                counts[assumption]["both"] += int(diag_rev and full_rev)
                if diag_rev and full_rev:
                    examples.append(ex)

    examples_sorted = sorted(examples, key=lambda e: e.full_fd_derivative)
    best = examples_sorted[:10]
    nontrivial_sorted = [
        ex for ex in examples_sorted if ex.assumption == "fixed_damping_per_inertia"
    ]
    best_nontrivial = nontrivial_sorted[:10]

    with (OUT / "inertia_reversal_counts.json").open("w", encoding="utf-8") as f:
        json.dump(counts, f, indent=2)
    with (OUT / "inertia_reversal_examples.json").open("w", encoding="utf-8") as f:
        json.dump([asdict(ex) for ex in best], f, indent=2)
    with (OUT / "inertia_reversal_nontrivial_examples.json").open("w", encoding="utf-8") as f:
        json.dump([asdict(ex) for ex in best_nontrivial], f, indent=2)
    with (OUT / "inertia_reversal_examples.csv").open("w", newline="", encoding="utf-8") as f:
        if best:
            writer = csv.DictWriter(f, fieldnames=list(asdict(best[0]).keys()))
            writer.writeheader()
            writer.writerows(asdict(ex) for ex in best)

    print(json.dumps(counts, indent=2))
    if best_nontrivial:
        ex = best_nontrivial[0]
        print("\nStrongest fixed_damping_per_inertia example:")
        print(json.dumps(asdict(ex), indent=2))
    else:
        print("\nNo fixed_damping_per_inertia case found with both diagonal and full-QEP reversal.")


if __name__ == "__main__":
    main()
