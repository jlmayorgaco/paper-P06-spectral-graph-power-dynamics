from __future__ import annotations

import numpy as np


def test_transverse_center_jordan_deflation():
    rng = np.random.default_rng(20260921)
    for _ in range(100):
        n_perp = int(rng.integers(2, 8))
        omega = float(rng.uniform(10.0, 100.0))
        b = rng.normal(size=(n_perp, n_perp))
        b -= (np.max(np.real(np.linalg.eigvals(b))) + 1.0) * np.eye(n_perp)
        s = rng.normal(size=(n_perp + 2, n_perp + 2))
        while abs(np.linalg.det(s)) < 1e-4:
            s = rng.normal(size=(n_perp + 2, n_perp + 2))
        j = np.zeros((n_perp + 2, n_perp + 2))
        j[0, 1] = omega
        j[2:, 2:] = b
        a = s @ j @ np.linalg.inv(s)
        rx = s[:, 0]
        w = s[:, 1]
        assert np.linalg.norm(a @ rx) / max(1.0, np.linalg.norm(a) * np.linalg.norm(rx)) < 1e-10
        assert np.linalg.norm(a @ w - omega * rx) / max(1.0, np.linalg.norm(a) * np.linalg.norm(w)) < 1e-10
        for z in (0.2 + 0.3j, 1.1 - 0.4j):
            lhs = np.linalg.det(z * np.eye(n_perp + 2) - a)
            rhs = z**2 * np.linalg.det(z * np.eye(n_perp) - b)
            assert abs(lhs - rhs) / max(1.0, abs(lhs)) < 1e-8

