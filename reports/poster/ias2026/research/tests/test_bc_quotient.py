"""BC01 synthetic obligations for T1/T1a and the four-state classification.

These are mathematical checks on constructed matrices, not benchmark results.
"""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.certification.classify import (
    BOUNDARY,
    STABLE,
    UNSTABLE,
    classify_spectrum,
)
from ibr_cycles.certification.symmetry import quotient, relocate


def jordan_example(eps: float) -> np.ndarray:
    """V1 of the theory note: symmetry e1, Jordan partner e2, physical eps."""

    return np.array([[0.0, 1.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, -2.0]])


def old_disc_count(a: np.ndarray, disc: float = 1e-3) -> int:
    v = np.linalg.eigvals(a)
    return int(np.count_nonzero((v.real > 0) & (np.abs(v) > disc)))


def test_quotient_of_jordan_example_is_diag():
    q = quotient(jordan_example(0.3), np.array([1.0, 0.0, 0.0]))
    assert q.symmetry_residual < 1e-15
    assert np.allclose(np.sort(np.linalg.eigvals(q.a_q).real), [-2.0, 0.3])


def test_tiny_physical_pole_stays_unstable():
    a = jordan_example(1e-7)
    q = quotient(a, np.array([1.0, 0.0, 0.0]))
    report = classify_spectrum(q.a_q, 0.0)
    assert report.status == UNSTABLE
    assert report.n_positive == 1
    # the excluded-disc rule of G1 would have hidden it
    assert old_disc_count(a) == 0


def test_exact_physical_zero_is_boundary_not_stable():
    q = quotient(jordan_example(0.0), np.array([1.0, 0.0, 0.0]))
    assert classify_spectrum(q.a_q, 0.0).status == BOUNDARY


def test_only_the_declared_symmetry_is_removed():
    """With one symmetry and a Jordan chain, the quotient keeps a physical zero."""

    a = jordan_example(0.0)
    q = quotient(a, np.array([1.0, 0.0, 0.0]))
    values = np.linalg.eigvals(q.a_q)
    assert np.count_nonzero(np.abs(values) < 1e-12) == 1
    assert classify_spectrum(q.a_q, 0.0).status == BOUNDARY


def _random_with_symmetry(seed: int, n: int = 7, r: int = 1):
    rng = np.random.default_rng(seed)
    basis, _ = np.linalg.qr(rng.standard_normal((n, n)))
    u = basis[:, :r]
    a_q = rng.standard_normal((n - r, n - r)) - 1.5 * np.eye(n - r)
    b_g = rng.standard_normal((r, n - r))
    t = np.block([[np.zeros((r, r)), b_g], [np.zeros((n - r, r)), a_q]])
    a = basis @ t @ basis.T
    # a non-orthonormal generator spanning the same subspace
    return a, u @ rng.standard_normal((r, r)), a_q


@pytest.mark.parametrize("beta", [0.1, 1.0, 10.0])
def test_relocation_moves_only_the_symmetry(beta):
    a, r, a_q = _random_with_symmetry(3)
    shifted = np.sort_complex(np.linalg.eigvals(relocate(a, r, beta)))
    expected = np.sort_complex(np.concatenate([[-beta], np.linalg.eigvals(a_q)]))
    assert np.allclose(shifted, expected, atol=1e-9)


def test_quotient_spectrum_is_basis_and_reference_independent():
    a, r, a_q = _random_with_symmetry(5)
    rng = np.random.default_rng(11)
    t = np.eye(a.shape[0]) + 0.3 * rng.standard_normal(a.shape)
    a2 = t @ a @ np.linalg.inv(t)
    r2 = t @ r
    q1 = np.sort_complex(np.linalg.eigvals(quotient(a, r).a_q))
    q2 = np.sort_complex(np.linalg.eigvals(quotient(a2, r2).a_q))
    assert np.allclose(q1, q2, atol=1e-8)
    assert np.allclose(q1, np.sort_complex(np.linalg.eigvals(a_q)), atol=1e-8)


def test_angle_reference_change_on_a_rotor_chain():
    """Absolute angles vs angles relative to machine 1: same physical spectrum."""

    rng = np.random.default_rng(2)
    n_m = 3
    # states (delta_1..3, omega_1..3); coupling depends on angle differences only
    k = rng.uniform(0.5, 2.0, size=(n_m, n_m))
    k = 0.5 * (k + k.T)
    np.fill_diagonal(k, 0.0)
    lap = np.diag(k.sum(1)) - k
    m = np.diag(rng.uniform(1, 3, n_m))
    d = np.diag(rng.uniform(0.1, 0.5, n_m))
    a = np.block(
        [
            [np.zeros((n_m, n_m)), np.eye(n_m)],
            [-np.linalg.solve(m, lap), -np.linalg.solve(m, d)],
        ]
    )
    r = np.concatenate([np.ones(n_m), np.zeros(n_m)])
    assert np.linalg.norm(a @ r) < 1e-12
    t = np.eye(2 * n_m)
    t[1:n_m, 0] = -1.0  # delta_i -> delta_i - delta_1 for i > 1
    a_rel = t @ a @ np.linalg.inv(t)
    q_abs = np.sort_complex(np.linalg.eigvals(quotient(a, r).a_q))
    q_rel = np.sort_complex(np.linalg.eigvals(quotient(a_rel, t @ r).a_q))
    assert np.allclose(q_abs, q_rel, atol=1e-10)
    assert classify_spectrum(quotient(a, r).a_q, 0.0).status == STABLE


def test_hidden_internal_mode_stays_in_the_full_ledger():
    """det P(s) = h(s) det T(s); a mode unobservable at the port is only in h."""

    a_d = np.diag([-1.0, 0.5])  # +0.5 is internal and unobservable
    b_d = np.array([[1.0], [1.0]])
    c_d = np.array([[1.0, 0.0]])
    d_d = np.array([[0.2]])
    y_net = np.array([[3.0]])

    def t_of(s):
        return y_net - d_d - c_d @ np.linalg.solve(s * np.eye(2) - a_d, b_d)

    def p_of(s):
        return np.block([[s * np.eye(2) - a_d, -b_d], [-c_d, y_net - d_d]])

    for s in (0.3 + 1j, 2.0, -0.4 + 0.2j):
        h = np.linalg.det(s * np.eye(2) - a_d)
        assert np.isclose(np.linalg.det(p_of(s)), h * np.linalg.det(t_of(s)))
    # the port determinant has neither a zero nor a pole at the hidden mode...
    assert 0.1 < abs(np.linalg.det(t_of(0.5 + 1e-9))) < 10.0
    # ...but the full model carries it: the eigenvalues of the closed loop
    closed = a_d + b_d @ np.linalg.solve(y_net - d_d, c_d)
    assert np.any(np.isclose(np.linalg.eigvals(closed), 0.5))
