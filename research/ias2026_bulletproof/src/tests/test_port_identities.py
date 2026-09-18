from __future__ import annotations

import numpy as np


def _blocks(rng: np.random.Generator, n: int):
    sizes = [int(rng.integers(1, 3)) for _ in range(n)]
    offsets = np.cumsum([0, *sizes])
    return sizes, offsets


def _matrix(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def test_determinant_and_local_collective_factorization_10000_cases():
    rng = np.random.default_rng(20260918)
    max_residual = 0.0
    for _ in range(10_000):
        n_blocks = int(rng.integers(1, 5))
        sizes, offsets = _blocks(rng, n_blocks)
        n = sum(sizes)
        d_blocks = [_matrix(rng, (s, s)) for s in sizes]
        d = np.zeros((n, n), dtype=complex)
        for i, (lo, hi) in enumerate(zip(offsets[:-1], offsets[1:])):
            d[lo:hi, lo:hi] = d_blocks[i]
        k = _matrix(rng, (n, n))
        # Keep the random cases away from accidental singular factors.
        k += 2.0 * np.eye(n)
        kd = np.zeros_like(k)
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            kd[lo:hi, lo:hi] = k[lo:hi, lo:hi]
        ko = k - kd
        l = np.eye(n, dtype=complex) + d @ kd
        q = np.linalg.solve(l, d @ ko)
        lhs = np.eye(n, dtype=complex) + d @ k
        rhs = l @ (np.eye(n, dtype=complex) + q)
        bridge_left = k @ np.linalg.inv(lhs)
        bridge_right = k @ np.linalg.inv(np.eye(n, dtype=complex) + q) @ np.linalg.inv(l)
        scale = max(1.0, np.linalg.norm(lhs), np.linalg.norm(rhs))
        max_residual = max(
            max_residual,
            np.linalg.norm(lhs - rhs) / scale,
            np.linalg.norm(bridge_left - bridge_right) / max(1.0, np.linalg.norm(bridge_left)),
        )
        det_residual = abs(np.linalg.det(lhs) - np.linalg.det(l) * np.linalg.det(np.eye(n) + q))
        max_residual = max(max_residual, det_residual / max(1.0, abs(np.linalg.det(lhs))))
    assert max_residual < 1e-10

