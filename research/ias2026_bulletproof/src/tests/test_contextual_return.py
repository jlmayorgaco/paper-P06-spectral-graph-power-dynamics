from __future__ import annotations

import numpy as np


def test_contextual_return_schur_identity_10000_cases():
    rng = np.random.default_rng(20260919)
    worst = 0.0
    for _ in range(10_000):
        ni = int(rng.integers(1, 4))
        nr = int(rng.integers(1, 7))
        # The theorem's Q is built from the off-diagonal interaction, so its
        # diagonal block is zero. Keeping this assumption explicit prevents a
        # false test of the stronger, generally invalid identity.
        qii = np.zeros((ni, ni), dtype=complex)
        qrr = rng.normal(size=(nr, nr)) + 1j * rng.normal(size=(nr, nr))
        qri = rng.normal(size=(nr, ni)) + 1j * rng.normal(size=(nr, ni))
        qir = rng.normal(size=(ni, nr)) + 1j * rng.normal(size=(ni, nr))
        qrr += 2.0 * np.eye(nr)
        q = np.block([[qii, qir], [qri, qrr]])
        a = np.eye(nr, dtype=complex) + qrr
        r = qir @ np.linalg.solve(a, qri)
        lhs = np.linalg.det(np.eye(ni + nr, dtype=complex) + q)
        rhs = np.linalg.det(a) * np.linalg.det(np.eye(ni, dtype=complex) - r)
        worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs)))
    assert worst < 1e-10
