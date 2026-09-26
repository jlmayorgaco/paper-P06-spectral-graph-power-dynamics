from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .green import ActionGreen

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class InteractionSplit:
    """Separation of a portfolio determinant into isolated and collective parts.

        det(I + M) = [product over a of det(I + M_aa)] * det(I + Q)

    The individual factor is exactly the product of the single-action
    determinant ratios, so the collective factor carries everything that no
    single action can explain. For blocks of rank above one the block-diagonal
    split is a choice, not a canonical object; the gauge invariant quantities
    are the cycle holonomies, not Q itself.
    """

    s: complex
    full: complex
    individual: complex
    collective: complex
    q: ComplexMatrix
    m: ComplexMatrix

    @property
    def factorization_error(self) -> float:
        scale = max(abs(self.full), 1e-300)
        return float(abs(self.full - self.individual * self.collective) / scale)

    @property
    def spectral_radius_q(self) -> float:
        return float(np.max(np.abs(np.linalg.eigvals(self.q))))

    @property
    def min_singular_value(self) -> float:
        identity = np.eye(self.q.shape[0], dtype=np.complex128)
        return float(np.linalg.svd(identity + self.q, compute_uv=False)[-1])

    @property
    def norm_q(self) -> float:
        return float(np.linalg.norm(self.q, 2))


def split_interaction(green: ActionGreen, s: complex) -> InteractionSplit:
    """Factor the portfolio determinant into individual and collective parts."""

    m = green.m(s)
    n = m.shape[0]
    identity = np.eye(n, dtype=np.complex128)
    total = identity + m
    d_self = np.zeros_like(total)
    individual = 1.0 + 0.0j
    for block in green.block_slices:
        diagonal = total[block, block]
        d_self[block, block] = diagonal
        individual *= complex(np.linalg.det(diagonal))
    q = np.linalg.solve(d_self, total) - identity
    return InteractionSplit(
        s=s,
        full=complex(np.linalg.det(total)),
        individual=individual,
        collective=complex(np.linalg.det(identity + q)),
        q=q,
        m=m,
    )
