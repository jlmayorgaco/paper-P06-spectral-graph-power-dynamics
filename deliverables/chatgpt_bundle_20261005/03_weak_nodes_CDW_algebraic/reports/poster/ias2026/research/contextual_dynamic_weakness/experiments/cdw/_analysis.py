"""Pure analysis helpers (no model evaluation): rankings, Shapley, submodularity, Gordan."""

from __future__ import annotations

import math
from itertools import combinations, permutations

import numpy as np


def optimal_linear_order(W: np.ndarray) -> tuple[list[int], float]:
    """Exact maximum of sum_{i before j} W[i, j] over permutations (DP over subsets)."""

    n = W.shape[0]
    full = (1 << n) - 1
    best = np.full(1 << n, -np.inf)
    arg = np.full(1 << n, -1, dtype=int)
    best[0] = 0.0
    # colsum[mask][j] = sum_{i in mask} W[i, j]
    for mask in range(1, full + 1):
        for j in range(n):
            if mask >> j & 1:
                prev = mask ^ (1 << j)
                gain = sum(W[i, j] for i in range(n) if prev >> i & 1)
                val = best[prev] + gain
                if val > best[mask]:
                    best[mask] = val
                    arg[mask] = j
    order = []
    mask = full
    while mask:
        j = arg[mask]
        order.append(j)
        mask ^= 1 << j
    order.reverse()
    return order, float(best[full])


def brute_linear_order(W):
    n = W.shape[0]
    bestv, besto = -np.inf, None
    for p in permutations(range(n)):
        v = sum(W[p[a], p[b]] for a in range(n) for b in range(a + 1, n))
        if v > bestv:
            bestv, besto = v, list(p)
    return besto, float(bestv)


def shapley(values: dict, players: tuple) -> dict:
    """Exact Shapley value of v(S) (keys: sorted tuples) for the given players."""

    n = len(players)
    out = {}
    for i in players:
        others = [p for p in players if p != i]
        tot = 0.0
        for r in range(n):
            w = math.factorial(r) * math.factorial(n - r - 1) / math.factorial(n)
            for s in combinations(others, r):
                s = tuple(sorted(s))
                si = tuple(sorted(s + (i,)))
                tot += w * (values[si] - values[s])
        out[i] = tot
    return out


def second_differences(values: dict, players: tuple):
    """d_ij(S) = v(S+i+j) - v(S+i) - v(S+j) + v(S) for all S, i<j not in S."""

    rows = []
    for r in range(len(players) - 1):
        for s in combinations(players, r):
            rest = [p for p in players if p not in s]
            for i, j in combinations(rest, 2):
                a = values.get(tuple(sorted(s)))
                b = values.get(tuple(sorted(s + (i,))))
                c = values.get(tuple(sorted(s + (j,))))
                d = values.get(tuple(sorted(s + (i, j))))
                if None in (a, b, c, d) or any(np.isnan(x) for x in (a, b, c, d)):
                    continue
                rows.append((tuple(sorted(s)), i, j, d - b - c + a))
    return rows


def is_submodular(values, players, tol=0.0) -> bool:
    """Direct definition check: Delta_i v(S) >= Delta_i v(T) for S subset T, i not in T."""

    for r in range(len(players) + 1):
        for T in combinations(players, r):
            T = tuple(sorted(T))
            for i in players:
                if i in T:
                    continue
                dT = values[tuple(sorted(T + (i,)))] - values[T]
                for q in range(len(T) + 1):
                    for S in combinations(T, q):
                        S = tuple(sorted(S))
                        dS = values[tuple(sorted(S + (i,)))] - values[S]
                        if dS < dT - tol:
                            return False
    return True


def _project_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection of v onto the probability simplex {lam >= 0, sum lam = 1}."""

    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    idx = np.arange(1, v.size + 1)
    rho = np.nonzero(u - css / idx > 0)[0][-1]
    theta = css[rho] / (rho + 1.0)
    return np.maximum(v - theta, 0.0)


def _gordan_projected_gradient(G: np.ndarray, iters: int = 20000):
    """Dependency-free fallback: FISTA (accelerated projected gradient) for min ||G^T lam||^2 on
    the simplex.

    Used only when every QP solver fails (numerical robustness, not a scientific choice: the
    objective, constraints and stopping rule are unchanged; see docs/CDW_PREREG_V1_DEVIATIONS.md).
    """

    m = G.shape[0]
    GGt = G @ G.T
    lip = max(2.0 * np.linalg.eigvalsh(GGt).max(), 1e-12)
    step = 1.0 / lip
    lam = np.full(m, 1.0 / m)
    y, t = lam.copy(), 1.0
    for _ in range(iters):
        grad = 2.0 * (GGt @ y)
        lam_new = _project_simplex(y - step * grad)
        t_new = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        y = lam_new + ((t - 1.0) / t_new) * (lam_new - lam)
        lam, t = lam_new, t_new
    return float(lam @ GGt @ lam), lam


def gordan(G: np.ndarray):
    """min || G^T lam ||^2 over the simplex (rows of G are gradients). Returns (value, lam).

    Tries CLARABEL, then SCS, then OSQP; falls back to a projected-gradient solve of the same
    convex problem if every cvxpy solver fails (ill-conditioned G). The problem, constraints and
    tolerance are unchanged across fallbacks.
    """

    import cvxpy as cp

    m = G.shape[0]
    if not np.all(np.isfinite(G)):
        return _gordan_projected_gradient(np.nan_to_num(G, nan=0.0, posinf=1e6, neginf=-1e6))
    lam = cp.Variable(m, nonneg=True)
    prob = cp.Problem(cp.Minimize(cp.sum_squares(G.T @ lam)), [cp.sum(lam) == 1])
    for solver in (cp.CLARABEL, cp.SCS, cp.OSQP):
        try:
            prob.solve(solver=solver)
            if prob.status in ("optimal", "optimal_inaccurate") and lam.value is not None:
                return float(prob.value), np.asarray(lam.value).ravel()
        except Exception:  # noqa: BLE001 - solver-specific failures, try the next one
            continue
    return _gordan_projected_gradient(G)


def sign_class(d, tau):
    if d <= -tau:
        return -1
    if d >= tau:
        return 1
    return 0


def entropy3(counts) -> float:
    tot = sum(counts)
    if tot == 0:
        return float("nan")
    return float(-sum((c / tot) * math.log2(c / tot) for c in counts if c > 0))


def kendall(a, b):
    from scipy.stats import kendalltau

    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return float(kendalltau(a[ok], b[ok]).statistic)


def spearman(a, b):
    from scipy.stats import spearmanr

    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
        return float("nan")
    return float(spearmanr(a[ok], b[ok]).statistic)


def topk_precision(pred, truth, k=5):
    """Overlap of the top-k by pred and by truth (both: larger = more)."""

    pred, truth = np.asarray(pred, float), np.asarray(truth, float)
    ok = np.isfinite(pred) & np.isfinite(truth)
    idx = np.where(ok)[0]
    if idx.size < k:
        return float("nan")
    tp = set(idx[np.argsort(-pred[idx])[:k]])
    tt = set(idx[np.argsort(-truth[idx])[:k]])
    return len(tp & tt) / k
