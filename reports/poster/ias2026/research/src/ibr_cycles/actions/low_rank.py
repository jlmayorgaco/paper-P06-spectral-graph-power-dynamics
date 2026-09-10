from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class ActionFactor:
    """Low-rank factorization of one action, written on the operator T(s).

    The convention throughout the package is

        T_S(s) = T_0(s) + sum over a in S of U_a C_a V_a^H
        T_0(s) = sI - A_0

    so for a state-matrix update dA_a the factored quantity is dT_a = -dA_a.
    The singular values, numerical rank and reconstruction error are stored
    rather than asserted: rank one is a measurement, not a declaration.
    """

    name: str
    u: ComplexMatrix
    c: ComplexMatrix
    v: ComplexMatrix
    singular_values: NDArray[np.float64]
    numerical_rank: int
    reconstruction_error: float
    tolerance: float

    @property
    def rank(self) -> int:
        return int(self.c.shape[0])

    def delta_t(self) -> ComplexMatrix:
        return self.u @ self.c @ self.v.conj().T

    def gauge_transform(self, basis: ComplexMatrix) -> ActionFactor:
        """Apply the internal basis change U -> U S and V^H -> S^-1 V^H.

        The factorization is not unique, so cycle holonomies must be invariant
        under this map. test_gauge_invariance.py checks exactly that.
        """

        inverse = np.linalg.inv(basis)
        return ActionFactor(
            name=self.name,
            u=self.u @ basis,
            c=inverse @ self.c @ basis,
            v=self.v @ inverse.conj().T,
            singular_values=self.singular_values,
            numerical_rank=self.numerical_rank,
            reconstruction_error=self.reconstruction_error,
            tolerance=self.tolerance,
        )


def factorize_delta(
    name: str,
    delta_t: NDArray[np.generic],
    *,
    relative_tolerance: float = 1e-10,
) -> ActionFactor:
    """SVD factorization of an operator update with a measured numerical rank."""

    matrix = np.asarray(delta_t, dtype=np.complex128)
    u, s, vh = np.linalg.svd(matrix)
    largest = float(s[0]) if s.size else 0.0
    cutoff = relative_tolerance * max(largest, 1e-300)
    rank = int(np.sum(s > cutoff))
    if rank == 0:
        raise ValueError(f"action {name!r} has a numerically zero update")
    u_r = u[:, :rank]
    v_r = vh[:rank, :].conj().T
    c_r = np.diag(s[:rank]).astype(np.complex128)
    reconstruction = u_r @ c_r @ v_r.conj().T
    error = float(np.linalg.norm(reconstruction - matrix) / max(largest, 1e-300))
    return ActionFactor(
        name=name,
        u=u_r,
        c=c_r,
        v=v_r,
        singular_values=np.asarray(s, dtype=np.float64),
        numerical_rank=rank,
        reconstruction_error=error,
        tolerance=relative_tolerance,
    )


def factorize_state_update(
    name: str,
    delta_a: NDArray[np.generic],
    *,
    relative_tolerance: float = 1e-10,
) -> ActionFactor:
    """Factorize a state-matrix action, converting dA into dT = -dA."""

    update = -np.asarray(delta_a, dtype=np.complex128)
    return factorize_delta(name, update, relative_tolerance=relative_tolerance)
