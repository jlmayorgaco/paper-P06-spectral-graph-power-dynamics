"""Stability verdicts that survive a system with no primary frequency control.

The IEEE-39 replacement model retains no governor, so the common angle and the
common frequency are undetermined and the reduced matrix carries a double
eigenvalue at the origin. Spectral abscissa is therefore not usable as a raw
metric: it is always about zero. Every verdict here excludes that reference pair
explicitly and asserts what it excluded, so a genuine mode drifting into the
exclusion band is reported rather than silently discarded.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..dynamics.modes import Mode, Spectrum

REFERENCE_THRESHOLD = 1e-3
OSCILLATORY_BAND_HZ = (0.1, 5.0)


@dataclass(frozen=True)
class Assessment:
    """Stability verdict of one linearized case, with its exclusion evidence."""

    spectral_abscissa: float
    critical_real: float
    critical_imag: float
    critical_frequency_hz: float
    critical_condition: float
    zeta_min: float
    zeta_min_frequency_hz: float
    zeta_min_real: float
    rhp_modes: int
    n_modes: int
    reference_modes: int
    reference_max_abs: float
    reference_participation: float
    oscillatory_modes: int

    @property
    def unstable(self) -> bool:
        return self.spectral_abscissa > 1e-6

    @property
    def reference_exclusion_is_clean(self) -> bool:
        """Exactly two excluded modes, dominated by machine angle and speed."""

        return self.reference_modes == 2 and self.reference_participation > 0.5


def _is_reference(mode: Mode, labels: tuple[str, ...]) -> float:
    weight = 0.0
    for value, name in zip(mode.participation, labels, strict=True):
        if name.startswith(("delta_sg", "omega_sg")):
            weight += float(value)
    return weight


def assess(
    spectrum: Spectrum,
    labels: tuple[str, ...],
    *,
    reference_threshold: float = REFERENCE_THRESHOLD,
    band_hz: tuple[float, float] = OSCILLATORY_BAND_HZ,
) -> Assessment:
    """Verdict over the non-reference modes of one spectrum."""

    reference = [m for m in spectrum.modes if abs(m.value) <= reference_threshold]
    dynamic = [m for m in spectrum.modes if abs(m.value) > reference_threshold]
    if not dynamic:
        raise ValueError("every mode fell inside the reference exclusion band")
    critical = max(dynamic, key=lambda m: m.real)
    band = [m for m in dynamic if band_hz[0] <= m.frequency_hz <= band_hz[1]]
    worst = min(band, key=lambda m: m.damping) if band else critical
    participation = (
        float(np.mean([_is_reference(m, labels) for m in reference]))
        if reference
        else 0.0
    )
    return Assessment(
        spectral_abscissa=critical.real,
        critical_real=critical.real,
        critical_imag=critical.imag,
        critical_frequency_hz=critical.frequency_hz,
        critical_condition=critical.condition,
        zeta_min=worst.damping,
        zeta_min_frequency_hz=worst.frequency_hz,
        zeta_min_real=worst.real,
        rhp_modes=sum(1 for m in dynamic if m.real > 1e-6),
        n_modes=len(spectrum.modes),
        reference_modes=len(reference),
        reference_max_abs=max((abs(m.value) for m in reference), default=0.0),
        reference_participation=participation,
        oscillatory_modes=len(band),
    )
