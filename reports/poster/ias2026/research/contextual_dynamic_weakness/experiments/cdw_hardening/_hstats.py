"""Statistics helpers (docs/CDW_HARDENING_STATISTICAL_PLAN.md): cluster bootstrap, sign test,
cluster sign-flip permutation, Holm. Intervals are resampling intervals over tested conditions."""

from __future__ import annotations

import numpy as np
from scipy.stats import binomtest

import _hinfra as HI

B = 10000


def cluster_boot(values, clusters=None, stat=np.median, seed=HI.SEED_BOOT, b=B):
    """Percentile bootstrap of stat over clusters (each resampled cluster contributes all its values).
    Returns (point, lo2.5, hi97.5, lo5)."""

    v = np.asarray(values, float)
    ok = np.isfinite(v)
    v = v[ok]
    if clusters is None:
        clusters = np.arange(v.size)
    else:
        clusters = np.asarray(clusters)[ok]
    ids = np.unique(clusters)
    groups = [v[clusters == c] for c in ids]
    rng = np.random.default_rng(seed)
    out = np.empty(b)
    for k in range(b):
        pick = rng.integers(len(groups), size=len(groups))
        out[k] = stat(np.concatenate([groups[j] for j in pick]))
    return float(stat(v)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)), float(np.percentile(out, 5))


def boot_diff_medians(a, b, seed=HI.SEED_BOOT, n=B):
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    rng = np.random.default_rng(seed)
    d = np.array([np.median(rng.choice(a, a.size)) - np.median(rng.choice(b, b.size)) for _ in range(n)])
    return float(np.median(a) - np.median(b)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def sign_test_less(values, thr):
    """k = #(values < thr) of n; exact one-sided binomial test against P(value < thr) <= 0.5."""

    v = np.asarray(values, float)
    v = v[np.isfinite(v) & (v != thr)]
    k = int((v < thr).sum())
    return k, int(v.size), float(binomtest(k, v.size, 0.5, alternative="greater").pvalue) if v.size else np.nan


def cluster_signflip_median(diffs, clusters, shift=0.0, seed=HI.SEED_PERM, n=B):
    """One-sided test that the median of (diffs - shift) > 0, flipping the sign of cluster means."""

    d = np.asarray(diffs, float) - shift
    c = np.asarray(clusters)
    ok = np.isfinite(d)
    d, c = d[ok], c[ok]
    ids = np.unique(c)
    means = np.array([d[c == k].mean() for k in ids])
    obs = np.median(means)
    rng = np.random.default_rng(seed)
    null = np.array([np.median(means * rng.choice([-1.0, 1.0], means.size)) for _ in range(n)])
    return float(obs), float((np.sum(null >= obs) + 1) / (n + 1))


def perm_spearman(x, y, seed=HI.SEED_PERM, n=B):
    from scipy.stats import spearmanr

    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 6 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan, int(x.size)
    rho = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y)).statistic for _ in range(n)])
    return rho, float((np.sum(np.abs(null) >= abs(rho)) + 1) / (n + 1)), int(x.size)


def holm(pvals: dict) -> dict:
    items = sorted(((k, v) for k, v in pvals.items() if np.isfinite(v)), key=lambda kv: kv[1])
    m = len(items)
    out, run = {}, 0.0
    for r, (k, p) in enumerate(items):
        run = max(run, min(1.0, (m - r) * p))
        out[k] = run
    for k, v in pvals.items():
        out.setdefault(k, np.nan)
    return out
