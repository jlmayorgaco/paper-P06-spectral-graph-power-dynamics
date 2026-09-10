from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..actions.green import ActionGreen

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class CycleHolonomy:
    """Closed intervention cycle a1 -> a2 -> ... -> ap -> a1.

        H_gamma = K_a1a2 K_a2a3 ... K_apa1

    Under an internal basis change H_gamma transforms by similarity, so trace,
    determinant and eigenvalues are gauge invariant and the entries are not.
    Only the invariants may be reported.
    """

    members: tuple[str, ...]
    s: complex
    matrix: ComplexMatrix

    @property
    def order(self) -> int:
        return len(self.members)

    @property
    def trace(self) -> complex:
        return complex(np.trace(self.matrix))

    @property
    def determinant(self) -> complex:
        return complex(np.linalg.det(self.matrix))

    @property
    def eigenvalues(self) -> NDArray[np.complex128]:
        return np.linalg.eigvals(self.matrix)

    @property
    def spectral_radius(self) -> float:
        return float(np.max(np.abs(self.eigenvalues)))

    @property
    def magnitude(self) -> float:
        """Scalar cycle score used for ranking; equals abs(q) for rank-one blocks."""

        return float(abs(self.trace)) if self.order > 0 else 0.0


def holonomy(green: ActionGreen, cycle: Sequence[str], s: complex) -> CycleHolonomy:
    """Evaluate the closed-cycle holonomy of an action cycle at frequency s."""

    members = tuple(cycle)
    if len(members) < 2:
        raise ValueError("a cycle needs at least two distinct actions")
    if len(set(members)) != len(members):
        raise ValueError("a simple cycle may not repeat an action")
    names = green.names
    index = {name: names.index(name) for name in members}
    k = green.k(s)
    product: ComplexMatrix | None = None
    for position, name in enumerate(members):
        following = members[(position + 1) % len(members)]
        factor = green.block(k, index[name], index[following])
        product = factor if product is None else product @ factor
    assert product is not None
    return CycleHolonomy(members=members, s=s, matrix=product)


def gauge_invariance_residual(
    green: ActionGreen,
    cycle: Sequence[str],
    s: complex,
    bases: Sequence[ComplexMatrix],
) -> dict[str, float]:
    """Change every internal basis and measure the drift of the invariants."""

    reference = holonomy(green, cycle, s)
    transformed = ActionGreen(
        a0=green.a0,
        factors=tuple(
            factor.gauge_transform(basis)
            for factor, basis in zip(green.factors, bases, strict=True)
        ),
    )
    moved = holonomy(transformed, cycle, s)
    trace_scale = max(abs(reference.trace), 1e-300)
    det_scale = max(abs(reference.determinant), 1e-300)
    eig_reference = np.sort_complex(reference.eigenvalues)
    eig_moved = np.sort_complex(moved.eigenvalues)
    eig_scale = max(float(np.abs(eig_reference).max()), 1e-300)
    return {
        "trace": float(abs(moved.trace - reference.trace) / trace_scale),
        "determinant": float(
            abs(moved.determinant - reference.determinant) / det_scale
        ),
        "eigenvalues": float(np.abs(eig_moved - eig_reference).max() / eig_scale),
    }
