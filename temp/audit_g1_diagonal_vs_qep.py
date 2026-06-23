from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def connected_laplacian(n: int, rng: np.random.Generator, extra_p: float = 0.12):
    W = np.zeros((n, n))
    # Random spanning tree first.
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
    zetas = [
        -lam.real / abs(lam)
        for lam in ev
        if abs(lam.imag) > 1e-7 and abs(lam) > 1e-9
    ]
    return min(zetas) if zetas else np.nan


def modal_data(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    Mh = np.diag(1 / np.sqrt(M))
    Lt = Mh @ L @ Mh
    Dt = np.diag(D / M)
    nu, Q = np.linalg.eigh(Lt)
    keep = nu > 1e-8
    nu = nu[keep]
    Q = Q[:, keep]
    delta = np.array([Q[:, k] @ Dt @ Q[:, k] for k in range(Q.shape[1])])
    zdiag = delta / (2 * np.sqrt(nu))
    return nu, Q, Dt, zdiag


def diag_margin(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    nu, _Q, _Dt, zdiag = modal_data(L, M, D)
    k = int(np.argmin(zdiag))
    return float(zdiag[k]), k, len(nu)


def reduced_qep_margin(L: np.ndarray, M: np.ndarray, D: np.ndarray, r: int = 6):
    n = len(M)
    Mh = np.diag(1 / np.sqrt(M))
    Lt = Mh @ L @ Mh
    Dt = np.diag(D / M)
    nu, Q = np.linalg.eigh(Lt)
    keep = np.where(nu > 1e-8)[0]
    zdiag = np.full(n, np.inf)
    for k in keep:
        zdiag[k] = (Q[:, k] @ Dt @ Q[:, k]) / (2 * np.sqrt(nu[k]))
    idx = np.argsort(zdiag)[: min(r, len(keep))]
    V = Q[:, idx]
    rr = V.shape[1]
    A = np.block(
        [
            [np.zeros((rr, rr)), np.eye(rr)],
            [-V.T @ Lt @ V, -V.T @ Dt @ V],
        ]
    )
    ev = np.linalg.eigvals(A)
    zetas = [
        -lam.real / abs(lam)
        for lam in ev
        if abs(lam.imag) > 1e-7 and abs(lam) > 1e-9
    ]
    return min(zetas) if zetas else np.nan


def rel_err(est: float, truth: float):
    return abs(est - truth) / abs(truth)


@dataclass
class Summary:
    name: str
    cases: int
    diag_med: float
    diag_p95: float
    qep_med: float
    qep_p95: float
    diag_better_frac: float
    critical_is_maxfreq_frac: float


def run_ensemble(name: str, make_case, cases: int = 300, r: int = 6):
    rng = np.random.default_rng(20260608)
    diag_errors = []
    qep_errors = []
    diag_better = 0
    maxfreq = 0
    valid = 0
    worst = None
    for _ in range(cases):
        L, M, D = make_case(rng)
        truth = exact_margin(L, M, D)
        if not np.isfinite(truth) or truth <= 1e-8:
            continue
        zd, k, m = diag_margin(L, M, D)
        zq = reduced_qep_margin(L, M, D, r=r)
        ed = rel_err(zd, truth)
        eq = rel_err(zq, truth)
        diag_errors.append(ed)
        qep_errors.append(eq)
        diag_better += ed <= eq
        maxfreq += k == (m - 1)
        valid += 1
        if worst is None or ed > worst[0]:
            worst = (ed, eq, truth, zd, zq, k, m)
    de = np.array(diag_errors)
    qe = np.array(qep_errors)
    return Summary(
        name=name,
        cases=valid,
        diag_med=float(np.median(de)),
        diag_p95=float(np.quantile(de, 0.95)),
        qep_med=float(np.median(qe)),
        qep_p95=float(np.quantile(qe, 0.95)),
        diag_better_frac=diag_better / valid,
        critical_is_maxfreq_frac=maxfreq / valid,
    ), worst


def main():
    def mild(rng):
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        d_over_m = np.exp(rng.uniform(math.log(0.08), math.log(0.8), n))
        D = M * d_over_m
        return L, M, D

    def weak_damping(rng):
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        d_over_m = np.exp(rng.uniform(math.log(0.006), math.log(0.08), n))
        D = M * d_over_m
        return L, M, D

    def proportional(rng):
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        D = 0.12 * M
        return L, M, D

    def strong_hetero(rng):
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(0.25), math.log(12.0), n))
        d_over_m = np.exp(rng.uniform(math.log(0.02), math.log(3.0), n))
        D = M * d_over_m
        return L, M, D

    def near_degenerate(rng):
        n = 24
        L = complete_laplacian(n, eps=0.015, rng=rng)
        M = np.ones(n)
        d_over_m = np.exp(rng.uniform(math.log(0.03), math.log(1.5), n))
        D = M * d_over_m
        return L, M, D

    for name, maker in [
        ("proportional_control", proportional),
        ("weak_damping_heterogeneity", weak_damping),
        ("random_mild_heterogeneity", mild),
        ("random_strong_heterogeneity", strong_hetero),
        ("near_degenerate_complete_graph", near_degenerate),
    ]:
        summary, worst = run_ensemble(name, maker)
        print(summary)
        print(
            "  worst_diag_case:",
            {
                "diag_err": worst[0],
                "qep_err": worst[1],
                "truth": worst[2],
                "diag": worst[3],
                "qep": worst[4],
                "diag_mode": worst[5],
                "num_modes": worst[6],
            },
        )

    # Exact repeated-frequency stress case.
    rng = np.random.default_rng(9)
    n = 18
    L = complete_laplacian(n, eps=0.0, rng=rng)
    M = np.ones(n)
    D = np.exp(rng.uniform(math.log(0.02), math.log(2.0), n))
    truth = exact_margin(L, M, D)
    zd, k, m = diag_margin(L, M, D)
    print(
        "exact_complete_graph_single_case",
        {
            "truth": truth,
            "diag": zd,
            "diag_err": rel_err(zd, truth),
            "qep_r6": reduced_qep_margin(L, M, D, r=6),
            "qep_r_all": reduced_qep_margin(L, M, D, r=n - 1),
            "diag_mode": k,
            "num_modes": m,
        },
    )


if __name__ == "__main__":
    main()
