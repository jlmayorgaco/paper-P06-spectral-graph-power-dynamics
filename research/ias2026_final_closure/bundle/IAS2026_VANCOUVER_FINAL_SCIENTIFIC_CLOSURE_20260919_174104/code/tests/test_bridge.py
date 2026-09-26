from __future__ import annotations

import numpy as np


def test_bridge_does_not_require_inverse_of_D():
    rng = np.random.default_rng(20260920)
    for _ in range(1000):
        n = int(rng.integers(2, 8))
        d = np.diag(np.r_[0.0, rng.normal(size=n - 1)]).astype(complex)
        k = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        k += 3.0 * np.eye(n)
        kd = np.diag(np.diag(k))
        l = np.eye(n, dtype=complex) + d @ kd
        q = np.linalg.solve(l, d @ (k - kd))
        lhs = k @ np.linalg.inv(np.eye(n) + d @ k)
        rhs = k @ np.linalg.inv(np.eye(n) + q) @ np.linalg.inv(l)
        assert np.linalg.norm(lhs - rhs) / max(1.0, np.linalg.norm(lhs)) < 1e-10

