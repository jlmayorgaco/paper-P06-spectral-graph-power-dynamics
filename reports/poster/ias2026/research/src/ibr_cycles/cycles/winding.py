from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Contour:
    """A closed, positively oriented contour sampled at ``points``."""

    name: str
    points: NDArray[np.complex128]
    metadata: dict[str, float]

    @property
    def samples(self) -> int:
        return int(self.points.size)

    def refined(self) -> Contour:
        raise NotImplementedError


@dataclass(frozen=True)
class Circle(Contour):
    def refined(self) -> Circle:
        return circle(
            complex(self.metadata["center_real"], self.metadata["center_imag"]),
            self.metadata["radius"],
            2 * self.samples,
        )


@dataclass(frozen=True)
class HalfPlane(Contour):
    def refined(self) -> HalfPlane:
        return right_half_plane(
            self.metadata["boundary"], self.metadata["radius"], 2 * self.samples
        )


def circle(center: complex, radius: float, samples: int = 2048) -> Circle:
    angles = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    return Circle(
        name="circle",
        points=center + radius * np.exp(1j * angles),
        metadata={
            "center_real": float(center.real),
            "center_imag": float(center.imag),
            "radius": float(radius),
        },
    )


def right_half_plane(
    boundary: float = -1e-3, radius: float = 100.0, samples: int = 8192
) -> HalfPlane:
    """Positively oriented D-contour enclosing the region Re(s) > boundary.

    Traversed counter-clockwise: down the arc and up the vertical segment, so a
    conjugate pair of enclosed zeros contributes +2, not +1. The count is per
    eigenvalue, never per conjugate pair.
    """

    arc_samples = samples // 2
    segment_samples = samples - arc_samples
    angles = np.linspace(-np.pi / 2.0, np.pi / 2.0, arc_samples, endpoint=False)
    arc = boundary + radius * np.exp(1j * angles)
    segment = boundary + 1j * np.linspace(
        radius, -radius, segment_samples, endpoint=False
    )
    return HalfPlane(
        name="right_half_plane",
        points=np.concatenate([arc, segment]),
        metadata={"boundary": float(boundary), "radius": float(radius)},
    )


@dataclass(frozen=True)
class Winding:
    """Argument-principle winding number with its numerical evidence."""

    contour: str
    samples: int
    raw: float
    minimum_magnitude: float
    refined_raw: float | None
    evaluations: int = 0
    max_phase_step: float = 0.0
    max_segment: float | None = None

    def resolved_for(self, clearance: float, *, factor: float = 4.0) -> bool:
        """Whether the sampling step is fine enough for the nearest zero or pole.

        A contour passing at distance ``clearance`` from a zero does almost all
        of its turning within that distance. If the step is not well below it,
        the excursion cancels between samples and the count is silently wrong.
        """

        return self.max_segment is not None and self.max_segment <= clearance / factor

    @property
    def nearest_integer(self) -> int:
        return int(round(self.raw))

    @property
    def integer_residual(self) -> float:
        return float(abs(self.raw - self.nearest_integer))

    @property
    def converged(self) -> bool:
        """Reliable only when it is integer AND stable under contour refinement."""

        if self.integer_residual > 1e-3:
            return False
        if self.refined_raw is None:
            return False
        return abs(self.refined_raw - self.raw) < 1e-3


def _segment_phase(
    function: Callable[[complex], complex],
    a: complex,
    b: complex,
    va: complex,
    vb: complex,
    threshold: float,
    depth: int,
    max_segment: float | None,
) -> tuple[float, float, int]:
    """Phase accumulated from a to b, bisecting until the step is resolved."""

    step = float(np.angle(vb / va))
    long_chord = max_segment is not None and abs(b - a) > max_segment
    if (abs(step) <= threshold and not long_chord) or depth == 0:
        return step, min(abs(va), abs(vb)), 0
    midpoint = 0.5 * (a + b)
    value = function(midpoint)
    if abs(value) == 0.0 or not np.isfinite(value):
        raise ValueError("function has a zero or a pole on the contour")
    left = _segment_phase(
        function, a, midpoint, va, value, threshold, depth - 1, max_segment
    )
    right = _segment_phase(
        function, midpoint, b, value, vb, threshold, depth - 1, max_segment
    )
    return (
        left[0] + right[0],
        min(left[1], right[1]),
        left[2] + right[2] + 1,
    )


def _accumulate(
    function: Callable[[complex], complex],
    contour: Contour,
    threshold: float,
    depth: int,
    max_segment: float | None = None,
) -> tuple[float, float, int]:
    values = np.array(
        [function(complex(z)) for z in contour.points], dtype=np.complex128
    )
    if float(np.abs(values).min()) <= 0.0 or not np.all(np.isfinite(values)):
        raise ValueError("function has a zero or a pole on the contour")
    total = 0.0
    minimum = float(np.abs(values).min())
    extra = 0
    n = values.size
    for index in range(n):
        following = (index + 1) % n
        phase, local_min, added = _segment_phase(
            function,
            complex(contour.points[index]),
            complex(contour.points[following]),
            complex(values[index]),
            complex(values[following]),
            threshold,
            depth,
            max_segment,
        )
        total += phase
        minimum = min(minimum, local_min)
        extra += added
    return total / (2.0 * np.pi), minimum, int(n + extra)


def winding_number(
    function: Callable[[complex], complex],
    contour: Contour,
    *,
    refine: bool = True,
    max_phase_step: float = np.pi / 4.0,
    max_depth: int = 26,
    max_segment: float | None = None,
) -> Winding:
    """Count zeros minus poles of ``function`` inside ``contour``.

    Sampling is adaptive: any segment whose phase turns by more than
    ``max_phase_step`` is bisected. Uniform sampling silently undercounts a
    lightly damped mode, because the determinant ratio does most of its turning
    within a distance of order ``zeta * abs(lambda)`` of the mode, which a fixed
    raster steps straight over. That failure returns a clean integer zero, so it
    cannot be caught by an integer-residual check.

    Phase-driven refinement alone is not sufficient either: when the contour
    passes within ``d`` of a zero and the sample spacing is much larger than
    ``d``, the excursion cancels between neighbouring samples and no bisection
    is triggered. Pass ``max_segment`` no larger than a quarter of the known
    clearance, and check ``Winding.resolved_for``.
    """

    raw, minimum, evaluations = _accumulate(
        function, contour, max_phase_step, max_depth, max_segment
    )
    refined = None
    if refine:
        refined = _accumulate(
            function, contour.refined(), max_phase_step, max_depth, max_segment
        )[0]
    return Winding(
        contour=contour.name,
        samples=contour.samples,
        raw=raw,
        minimum_magnitude=minimum,
        refined_raw=refined,
        evaluations=evaluations,
        max_phase_step=max_phase_step,
        max_segment=max_segment,
    )


@dataclass(frozen=True)
class Admissibility:
    """Whether a region Omega supports a provenance statement at all.

    The individual/collective winding split is only interpretable when Omega
    contains the newly created mode and excludes every base mode and every mode
    of a lower-order portfolio. Otherwise both windings are still arithmetically
    correct and the interpretation silently changes: a contour that also
    encloses a single-action mode moves the count from the collective factor to
    the individual one. This test makes that failure loud.
    """

    admissible: bool
    enclosed_new: int
    enclosed_base: int
    enclosed_lower_order: int
    clearance: float
    reason: str


def _inside(contour: Contour, values: NDArray[np.complex128]) -> NDArray[np.bool_]:
    path = contour.points
    closed = np.concatenate([path, path[:1]])
    counts = []
    for value in values:
        offsets = closed - value
        turning = np.sum(np.angle(offsets[1:] / offsets[:-1]))
        counts.append(abs(turning) > np.pi)
    return np.array(counts, dtype=bool)


def check_admissibility(
    contour: Contour,
    *,
    portfolio_modes: NDArray[np.complex128],
    base_modes: NDArray[np.complex128],
    lower_order_modes: NDArray[np.complex128],
) -> Admissibility:
    """Verify that Omega isolates the created mode from all reference modes."""

    new_inside = int(np.sum(_inside(contour, portfolio_modes)))
    base_inside = int(np.sum(_inside(contour, base_modes)))
    lower_inside = int(np.sum(_inside(contour, lower_order_modes)))
    reference = np.concatenate([base_modes, lower_order_modes])
    distances = np.abs(contour.points[:, None] - reference[None, :])
    clearance = float(distances.min()) if reference.size else float("inf")
    if new_inside == 0:
        return Admissibility(
            False,
            new_inside,
            base_inside,
            lower_inside,
            clearance,
            "Omega contains no portfolio mode",
        )
    if base_inside or lower_inside:
        return Admissibility(
            False,
            new_inside,
            base_inside,
            lower_inside,
            clearance,
            "Omega also encloses base or lower-order modes, so the "
            "individual/collective split is not interpretable",
        )
    return Admissibility(
        True, new_inside, base_inside, lower_inside, clearance, "admissible"
    )


def spectra_union(spectra: Sequence[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    if not spectra:
        return np.zeros(0, dtype=np.complex128)
    return np.concatenate([np.asarray(s, dtype=np.complex128) for s in spectra])
