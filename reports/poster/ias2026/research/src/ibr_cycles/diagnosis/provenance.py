from __future__ import annotations

from dataclasses import dataclass

from ..actions.green import ActionGreen
from ..actions.interaction import split_interaction
from ..cycles.winding import Admissibility, Contour, Winding, winding_number


@dataclass(frozen=True)
class Provenance:
    """Attribution of the modes inside Omega to isolated or collective effects.

    The three windings satisfy full = individual + collective exactly, because
    the determinant factorizes. The scientific content is which term carries the
    count, and that statement is only meaningful on an admissible Omega.
    """

    full: Winding
    individual: Winding
    collective: Winding
    admissibility: Admissibility

    @property
    def consistent(self) -> bool:
        return abs(self.full.raw - self.individual.raw - self.collective.raw) < 1e-6

    @property
    def converged(self) -> bool:
        return (
            self.full.converged
            and self.individual.converged
            and self.collective.converged
        )

    @property
    def verdict(self) -> str:
        if not self.admissibility.admissible:
            return "INADMISSIBLE_REGION"
        if not self.converged:
            return "NUMERICALLY_UNRESOLVED"
        if not self.consistent:
            return "FACTORIZATION_INCONSISTENT"
        if (
            self.collective.nearest_integer != 0
            and self.individual.nearest_integer == 0
        ):
            return "COLLECTIVE"
        if (
            self.individual.nearest_integer != 0
            and self.collective.nearest_integer == 0
        ):
            return "INDIVIDUAL"
        return "MIXED"

    def as_row(self) -> dict[str, float | int | str]:
        return {
            "verdict": self.verdict,
            "winding_full": self.full.raw,
            "winding_individual": self.individual.raw,
            "winding_collective": self.collective.raw,
            "integer_full": self.full.nearest_integer,
            "integer_individual": self.individual.nearest_integer,
            "integer_collective": self.collective.nearest_integer,
            "residual_full": self.full.integer_residual,
            "min_magnitude_full": self.full.minimum_magnitude,
            "admissible": int(self.admissibility.admissible),
            "clearance": self.admissibility.clearance,
            "admissibility_reason": self.admissibility.reason,
        }


def decompose_winding(
    green: ActionGreen, contour: Contour, admissibility: Admissibility
) -> Provenance:
    """Split the portfolio winding into individual and collective contributions."""

    def full(s: complex) -> complex:
        return split_interaction(green, s).full

    def individual(s: complex) -> complex:
        return split_interaction(green, s).individual

    def collective(s: complex) -> complex:
        return split_interaction(green, s).collective

    return Provenance(
        full=winding_number(full, contour),
        individual=winding_number(individual, contour),
        collective=winding_number(collective, contour),
        admissibility=admissibility,
    )
