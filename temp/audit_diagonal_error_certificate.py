from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def connected_laplacian(n: int, rng: np.random.Generator, extra_p: float = 0.10):
    W = np.zeros((n, n))
    for i in range(1, n):
        j = int(rng.integers(0, i))
        w = float(np.exp(rng.uniform(math.log(0.4), math.log(4.0))))
        W[i, j] = W[j, i] = w
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] == 0 and rng.random() < extra_p:
                w = float(np.exp(rng.uniform(math.log(0.2), math.log(6.0))))
                W[i, j] = W[j, i] = w
    return np.diag(W.sum(axis=1)) - W


def complete_laplacian(n: int, eps: float, rng: np.random.Generator):
    W = np.ones((n, n)) - np.eye(n)
    if eps:
        noise = rng.normal(0, eps, size=(n, n))
        noise = (noise + noise.T) / 2
        W = np.maximum(0.05, W * (1 + noise))
        np.fill_diagonal(W, 0)
    return np.diag(W.sum(axis=1)) - W


def exact_margin(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    n = len(M)
    A = np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-np.diag(1 / M) @ L, -np.diag(D / M)],
        ]
    )
    ev = np.linalg.eigvals(A)
    zetas = np.array(
        [
            -lam.real / abs(lam)
            for lam in ev
            if abs(lam.imag) > 1e-7 and abs(lam) > 1e-9
        ]
    )
    return float(np.min(zetas)) if len(zetas) else np.nan


def modal_objects(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    Mh = np.diag(1 / np.sqrt(M))
    Lt = Mh @ L @ Mh
    Dt = np.diag(D / M)
    nu_all, Q_all = np.linalg.eigh(Lt)
    keep = nu_all > 1e-8
    nu = nu_all[keep]
    Q = Q_all[:, keep]
    Gamma = Q.T @ Dt @ Q
    delta = np.diag(Gamma).copy()
    off = Gamma - np.diag(delta)
    zdiag = delta / (2 * np.sqrt(nu))
    return nu, Gamma, off, zdiag


def diagonal_poles(nu: np.ndarray, delta: np.ndarray):
    poles = []
    mode_ids = []
    for k, (vk, dk) in enumerate(zip(nu, delta)):
        disc = vk - 0.25 * dk * dk
        if disc <= 1e-12:
            continue
        wd = math.sqrt(disc)
        poles.append(complex(-0.5 * dk, wd))
        poles.append(complex(-0.5 * dk, -wd))
        mode_ids.extend([k, k])
    return np.array(poles, dtype=complex), np.array(mode_ids, dtype=int)


def diag_margin(nu: np.ndarray, zdiag: np.ndarray):
    under = zdiag < 1.0
    if not np.any(under):
        return np.nan, -1
    idxs = np.where(under)[0]
    k = int(idxs[np.argmin(zdiag[under])])
    return float(zdiag[k]), k


def bauer_fike_bound(nu: np.ndarray, Gamma: np.ndarray, off: np.ndarray):
    delta = np.diag(Gamma)
    poles, mode_ids = diagonal_poles(nu, delta)
    if len(poles) == 0:
        return math.inf, math.inf

    m = len(nu)
    cols = []
    for s, k in zip(poles, mode_ids):
        v = np.zeros(m)
        v[k] = 1.0
        cols.append(np.r_[v, s * v])
    V = np.column_stack(cols)
    if V.shape[0] != V.shape[1]:
        return math.inf, math.inf
    try:
        kappa = float(np.linalg.cond(V))
    except np.linalg.LinAlgError:
        return math.inf, math.inf
    rho = kappa * float(np.linalg.norm(off, 2))
    min_abs = float(np.min(np.abs(poles)))
    if not np.isfinite(rho) or rho >= min_abs:
        return rho, math.inf
    return rho, rho / (min_abs - rho)


def row_gershgorin_radius(abs_s: float, sep: float, row_sum: float):
    # Smallest positive r satisfying r(sep-r) >= (|s|+r) row_sum.
    # This is a sufficient local diagonal-dominance radius. If the
    # discriminant is negative, the bound is not informative.
    if row_sum <= 1e-14:
        return 0.0
    a = 1.0
    b = -(sep - row_sum)
    c = abs_s * row_sum
    disc = b * b - 4 * a * c
    if disc <= 0:
        return math.inf
    r = (-b - math.sqrt(disc)) / (2 * a)
    return float(r) if r >= 0 else math.inf


def row_certificate_bound(nu: np.ndarray, Gamma: np.ndarray, off: np.ndarray):
    delta = np.diag(Gamma)
    poles, mode_ids = diagonal_poles(nu, delta)
    if len(poles) == 0:
        return math.inf, math.inf, False

    radii = []
    for idx, (s, k) in enumerate(zip(poles, mode_ids)):
        other = np.delete(poles, idx)
        sep = float(np.min(np.abs(s - other))) if len(other) else math.inf
        row = float(np.sum(np.abs(off[k, :])))
        r = row_gershgorin_radius(abs(s), sep, row)
        radii.append(r)

    radii = np.array(radii)
    if not np.all(np.isfinite(radii)):
        return math.inf, math.inf, False

    # Conservative isolation check: if local disks overlap, use the bound only
    # as a trigger, not as a certified mode-by-mode error bound.
    isolated = True
    for i in range(len(poles)):
        for j in range(i + 1, len(poles)):
            if abs(poles[i] - poles[j]) <= radii[i] + radii[j]:
                isolated = False
                break
        if not isolated:
            break

    min_abs_after = float(np.min(np.abs(poles) - radii))
    if min_abs_after <= 0:
        return float(np.max(radii)), math.inf, False
    zeta_bound = float(np.max(radii / (np.abs(poles) - radii)))
    return float(np.max(radii)), zeta_bound, isolated


def second_order_qep_indicator(nu: np.ndarray, Gamma: np.ndarray, off: np.ndarray):
    """Second-order pole-shift indicator for off-diagonal modal damping.

    In the diagonal modal QEP, the first-order shift from E=off(Gamma) is zero
    because E_kk=0. Eliminating the other modal coordinates gives the
    second-order estimate

        |Delta s_k| <= |s_k|^2 / |2s_k+delta_k|
                      sum_l |E_kl|^2 / |p_l(s_k)|.

    This is an a-priori indicator, not a rigorous global certificate unless
    supplemented by a small-gain condition. It is included because it has the
    right mechanism: squared coupling and inverse modal separation.
    """
    delta = np.diag(Gamma)
    zdiag = delta / (2 * np.sqrt(nu))
    mode_bounds = []
    for k, (vk, dk, zk) in enumerate(zip(nu, delta, zdiag)):
        disc = vk - 0.25 * dk * dk
        if disc <= 1e-12 or zk >= 1.0:
            continue
        s = complex(-0.5 * dk, math.sqrt(disc))
        denom = abs(2 * s + dk)
        if denom <= 1e-14:
            mode_bounds.append(math.inf)
            continue
        acc = 0.0
        for ell in range(len(nu)):
            if ell == k:
                continue
            p_ell = s * s + delta[ell] * s + nu[ell]
            if abs(p_ell) <= 1e-14:
                acc = math.inf
                break
            acc += (abs(off[k, ell]) ** 2) / abs(p_ell)
        pole_radius = (abs(s) ** 2 / denom) * acc
        if not np.isfinite(pole_radius) or pole_radius >= abs(s):
            zeta_bound = math.inf
        else:
            zeta_bound = pole_radius / (abs(s) - pole_radius)
        mode_bounds.append(float(zeta_bound))
    if not mode_bounds:
        return math.inf
    return float(np.max(mode_bounds))


def coupling_metrics(nu: np.ndarray, off: np.ndarray):
    roots = np.sqrt(nu)
    gaps = np.abs(roots[:, None] - roots[None, :])
    mask = ~np.eye(len(nu), dtype=bool)
    min_gap = float(np.min(gaps[mask])) if len(nu) > 1 else math.inf
    eta = float(np.max(np.abs(off)[mask] / (gaps[mask] + 1e-12))) if len(nu) > 1 else 0.0
    return {
        "off_norm_2": float(np.linalg.norm(off, 2)),
        "off_norm_fro": float(np.linalg.norm(off, "fro")),
        "min_sqrt_gap": min_gap,
        "eta_max": eta,
        "off_over_gap": float(np.linalg.norm(off, 2) / (min_gap + 1e-12)),
    }


@dataclass
class CaseResult:
    family: str
    truth: float
    diag: float
    abs_error: float
    rel_error: float
    bf_pole_radius: float
    bf_zeta_bound: float
    row_pole_radius: float
    row_zeta_bound: float
    row_isolated: bool
    second_order_zeta_bound: float
    off_norm_2: float
    off_norm_fro: float
    min_sqrt_gap: float
    eta_max: float
    off_over_gap: float


def evaluate_case(family: str, L: np.ndarray, M: np.ndarray, D: np.ndarray):
    truth = exact_margin(L, M, D)
    if not np.isfinite(truth) or truth <= 1e-10:
        return None
    nu, Gamma, off, zdiag = modal_objects(L, M, D)
    diag, _ = diag_margin(nu, zdiag)
    if not np.isfinite(diag):
        return None
    bf_r, bf_z = bauer_fike_bound(nu, Gamma, off)
    row_r, row_z, row_iso = row_certificate_bound(nu, Gamma, off)
    second_order_zeta = second_order_qep_indicator(nu, Gamma, off)
    cm = coupling_metrics(nu, off)
    err = abs(diag - truth)
    return CaseResult(
        family=family,
        truth=truth,
        diag=diag,
        abs_error=float(err),
        rel_error=float(err / abs(truth)),
        bf_pole_radius=bf_r,
        bf_zeta_bound=bf_z,
        row_pole_radius=row_r,
        row_zeta_bound=row_z,
        row_isolated=row_iso,
        second_order_zeta_bound=second_order_zeta,
        **cm,
    )


def make_family(name: str, rng: np.random.Generator):
    if name == "proportional":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        D = 0.12 * M
        return L, M, D
    if name == "weak_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        ratio = np.exp(rng.uniform(math.log(0.006), math.log(0.08), n))
        D = M * ratio
        return L, M, D
    if name == "mild_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        ratio = np.exp(rng.uniform(math.log(0.08), math.log(0.8), n))
        D = M * ratio
        return L, M, D
    if name == "strong_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(0.25), math.log(12.0), n))
        ratio = np.exp(rng.uniform(math.log(0.02), math.log(3.0), n))
        D = M * ratio
        return L, M, D
    if name == "near_degenerate":
        n = 24
        L = complete_laplacian(n, eps=0.015, rng=rng)
        M = np.ones(n)
        ratio = np.exp(rng.uniform(math.log(0.03), math.log(1.5), n))
        D = M * ratio
        return L, M, D
    if name.startswith("rho_"):
        rho = float(name.split("_", 1)[1])
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        base = 0.12
        hetero = np.exp(rng.uniform(math.log(0.01), math.log(1.2), n))
        ratio = (1 - rho) * base + rho * hetero
        D = M * ratio
        return L, M, D
    raise ValueError(name)


def summarize(results: list[CaseResult]):
    arr = {k: np.array([getattr(r, k) for r in results], dtype=float) for k in [
        "abs_error",
        "rel_error",
        "bf_zeta_bound",
        "row_zeta_bound",
        "second_order_zeta_bound",
        "eta_max",
        "off_over_gap",
    ]}
    finite_bf = np.isfinite(arr["bf_zeta_bound"])
    finite_row = np.isfinite(arr["row_zeta_bound"])
    finite_so = np.isfinite(arr["second_order_zeta_bound"])
    isolated = np.array([r.row_isolated for r in results], dtype=bool)

    def q(x, p):
        return float(np.quantile(x, p)) if len(x) else math.nan

    def coverage(bound):
        finite = np.isfinite(bound)
        if not np.any(finite):
            return math.nan
        return float(np.mean(bound[finite] >= arr["abs_error"][finite]))

    def tightness(bound):
        finite = np.isfinite(bound) & (arr["abs_error"] > 1e-12)
        if not np.any(finite):
            return math.inf
        return float(np.median(bound[finite] / arr["abs_error"][finite]))

    return {
        "cases": len(results),
        "abs_error_median": q(arr["abs_error"], 0.5),
        "abs_error_p95": q(arr["abs_error"], 0.95),
        "rel_error_median": q(arr["rel_error"], 0.5),
        "rel_error_p95": q(arr["rel_error"], 0.95),
        "bf_finite_fraction": float(np.mean(finite_bf)),
        "bf_coverage_when_finite": coverage(arr["bf_zeta_bound"]),
        "bf_tightness_median": tightness(arr["bf_zeta_bound"]),
        "row_finite_fraction": float(np.mean(finite_row)),
        "row_isolated_fraction": float(np.mean(isolated)),
        "row_coverage_when_finite": coverage(arr["row_zeta_bound"]),
        "row_tightness_median": tightness(arr["row_zeta_bound"]),
        "second_order_finite_fraction": float(np.mean(finite_so)),
        "second_order_coverage_when_finite": coverage(arr["second_order_zeta_bound"]),
        "second_order_tightness_median": tightness(arr["second_order_zeta_bound"]),
        "second_order_bound_median": q(arr["second_order_zeta_bound"][finite_so], 0.5) if np.any(finite_so) else math.inf,
        "second_order_bound_p95": q(arr["second_order_zeta_bound"][finite_so], 0.95) if np.any(finite_so) else math.inf,
        "eta_median": q(arr["eta_max"], 0.5),
        "eta_p95": q(arr["eta_max"], 0.95),
        "off_over_gap_median": q(arr["off_over_gap"], 0.5),
        "off_over_gap_p95": q(arr["off_over_gap"], 0.95),
    }


def main():
    rng = np.random.default_rng(20260608)
    families = [
        "proportional",
        "weak_ibr",
        "mild_ibr",
        "strong_ibr",
        "near_degenerate",
        "rho_0.00",
        "rho_0.20",
        "rho_0.50",
        "rho_0.80",
    ]
    all_results = []
    report = {}
    for family in families:
        results = []
        for _ in range(250):
            L, M, D = make_family(family, rng)
            r = evaluate_case(family, L, M, D)
            if r is not None:
                results.append(r)
        all_results.extend(results)
        report[family] = summarize(results)

    out = {
        "interpretation": {
            "abs_error": "|zeta_diag - zeta_full|",
            "bf_zeta_bound": "global Bauer-Fike companion bound translated through zeta(s)",
            "row_zeta_bound": "local row/Gershgorin-style bound; useful only when finite and isolated",
            "row_isolated": "local modal disks do not overlap; otherwise the row bound is only a trigger",
            "second_order_zeta_bound": "QEP-specific second-order indicator; not yet a rigorous certificate",
        },
        "families": report,
        "sample_results": [asdict(r) for r in all_results[:20]],
    }
    out_path = Path("temp") / "diagonal_error_certificate_report.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
