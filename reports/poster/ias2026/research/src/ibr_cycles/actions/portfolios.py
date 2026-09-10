from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from itertools import combinations


def subsets(
    names: Sequence[str], *, include_empty: bool = True
) -> Iterator[tuple[str, ...]]:
    """Every subset of an action library, in nondecreasing size order."""

    start = 0 if include_empty else 1
    for size in range(start, len(names) + 1):
        yield from combinations(names, size)


@dataclass(frozen=True)
class PortfolioOutcome:
    """Stability verdict of one portfolio at one operating point."""

    members: tuple[str, ...]
    spectral_abscissa: float
    critical_real: float
    critical_imag: float
    damping: float
    frequency_hz: float
    rhp_modes: int
    eigenvalue_condition: float

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def stable(self) -> bool:
        return self.spectral_abscissa < 0.0


def minimum_destabilizing_order(outcomes: Sequence[PortfolioOutcome]) -> int | None:
    """kappa_Omega: the smallest portfolio size that is unstable.

    kappa is a property of the frozen action amplitudes, not of the system:
    scaling every action changes it. Callers must report it together with the
    amplitude vector that produced it.
    """

    unstable = [outcome.size for outcome in outcomes if not outcome.stable]
    return min(unstable) if unstable else None


def lower_order_safe(outcomes: Sequence[PortfolioOutcome], size: int) -> bool:
    """True when every portfolio strictly smaller than ``size`` is stable."""

    return all(outcome.stable for outcome in outcomes if outcome.size < size)
