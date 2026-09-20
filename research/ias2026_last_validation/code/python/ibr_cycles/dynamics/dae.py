from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

RealVector = NDArray[np.float64]
RealMatrix = NDArray[np.float64]


@runtime_checkable
class SemiExplicitDAE(Protocol):
    """A nonlinear semi-explicit index-1 DAE.

    ``xdot = f(x, z, theta)`` and ``0 = g(x, z, theta)``.
    """

    n_x: int
    n_z: int

    def f(
        self, x: RealVector, z: RealVector, theta: dict[str, float]
    ) -> RealVector: ...

    def g(
        self, x: RealVector, z: RealVector, theta: dict[str, float]
    ) -> RealVector: ...


@dataclass(frozen=True)
class Residual:
    """Stacked equilibrium residual and its component norms."""

    f: RealVector
    g: RealVector

    @property
    def stacked(self) -> RealVector:
        return np.concatenate([self.f, self.g])

    @property
    def norm_f(self) -> float:
        return float(np.linalg.norm(self.f, np.inf))

    @property
    def norm_g(self) -> float:
        return float(np.linalg.norm(self.g, np.inf))

    @property
    def norm(self) -> float:
        return max(self.norm_f, self.norm_g)


def residual(
    model: SemiExplicitDAE,
    x: RealVector,
    z: RealVector,
    theta: dict[str, float],
) -> Residual:
    return Residual(
        f=np.asarray(model.f(x, z, theta), dtype=np.float64),
        g=np.asarray(model.g(x, z, theta), dtype=np.float64),
    )


@dataclass(frozen=True)
class LinearSystem:
    """Reduced dynamic matrix ``Ared`` with state labels.

    ``Ared`` is the local small-signal ground truth after index-1 elimination of
    the algebraic variables.
    """

    A: RealMatrix
    labels: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.A.shape[0] != self.A.shape[1]:
            raise ValueError("Ared must be square")
        if len(self.labels) != self.A.shape[0]:
            raise ValueError("one label per differential state is required")

    @property
    def n(self) -> int:
        return int(self.A.shape[0])

    def index(self, label: str) -> int:
        return self.labels.index(label)

    def indices(self, labels: tuple[str, ...]) -> list[int]:
        return [self.index(name) for name in labels]
