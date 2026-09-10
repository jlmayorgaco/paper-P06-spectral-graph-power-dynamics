"""Frozen v2C definitions, shared so that every v2C script uses the same code.

The single methodological change from v2B lives here: the stability endpoint is
the envelope of the inter-area modal FAMILY, not one selected descendant. The
anchor rule, the band, the MAC floor and the core are inherited unchanged from
the frozen v2B lineage, so exactly one thing differs between the campaigns.

Protocol: configs/ias2026/v2c_modal_family_protocol.yaml
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ibr_cycles.dynamics.modal_family import (
    ModalFamily,
    band_candidates,
    descendant_family,
)
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import Spectrum, eigen_analysis
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

CORE = (30, 33, 35, 37)
BAND_HZ = (0.3, 1.5)
"""The frozen inter-area band, unchanged since trackA_final_v1."""

FAMILY_THRESHOLD = 0.80
"""Shape-overlap floor for family membership. The same 0.80 MAC floor v2B used."""

ZERO_MODE = 1e-3


def machine_labels(labels) -> set[str]:
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def nominal_reference() -> tuple[dict[str, complex], set[str]]:
    """The frozen inter-area shape, from the flagship at the nominal point.

    Inherited verbatim from v2B. The anchor rule was never the defective part;
    the defect was selecting a single descendant of the anchor.
    """

    case = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    mode = max(
        (m for m in eigen_analysis(case.system.A).modes if abs(m.value) > ZERO_MODE),
        key=lambda m: m.real,
    )
    shape = {
        name: mode.right[i]
        for i, name in enumerate(case.system.labels)
        if name.startswith(("delta_sg", "omega_sg"))
    }
    return shape, machine_labels(case.system.labels)


def base_anchor(spectrum: Spectrum, case, nominal):
    """The base-case inter-area mode: the band mode closest in shape to nominal."""

    shape, names = nominal
    common = sorted(set(names) & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    reference = np.array([shape[n] for n in common])
    index = {n: i for i, n in enumerate(case.system.labels)}
    positions = [index[n] for n in common]
    candidates = band_candidates(spectrum, BAND_HZ)
    if not candidates:
        return None
    return max(
        candidates, key=lambda m: modal_assurance(reference, m.right[positions])
    )


@dataclass(frozen=True)
class FamilyReading:
    """One case read through the frozen family observable."""

    family: ModalFamily
    argmax_real: float
    argmax_frequency_hz: float
    argmax_overlap: float

    @property
    def alpha(self) -> float:
        return self.family.alpha

    @property
    def frequency_worst_hz(self) -> float:
        return self.family.frequency_worst_hz


def read_family(
    anchor, base_case, case, spectrum, *, common_labels=None
) -> FamilyReading | None:
    """Family envelope of ``case``, against the frozen anchor of ``base_case``.

    Returns ``None`` for a tracking failure. The caller must classify that as
    TRACKING_FAILURE and never fall back to the nearest frequency.

    ``common_labels`` pins the comparison coordinates. Over a subset lattice the
    survivors differ from subset to subset, so the overlap threshold would mean
    different things at different nodes unless one fixed basis is used for all of
    them. The flagship survivors are contained in every proper subset's
    survivors, which makes them the natural choice.
    """

    common = sorted(
        machine_labels(base_case.system.labels) & machine_labels(case.system.labels)
        if common_labels is None
        else set(common_labels) & machine_labels(case.system.labels)
    )
    if len(common) < 4:
        return None
    base_index = {n: i for i, n in enumerate(base_case.system.labels)}
    reference = anchor.right[[base_index[n] for n in common]]
    norm = float(np.linalg.norm(reference))
    if norm == 0.0:
        return None
    reference = reference / norm
    index = {n: i for i, n in enumerate(case.system.labels)}
    positions = [index[n] for n in common]
    family = descendant_family(
        reference,
        positions,
        spectrum,
        band_hz=BAND_HZ,
        threshold=FAMILY_THRESHOLD,
    )
    if family is None:
        return None
    candidates = band_candidates(spectrum, BAND_HZ)
    best = max(candidates, key=lambda m: modal_assurance(reference, m.right[positions]))
    return FamilyReading(
        family=family,
        argmax_real=best.real,
        argmax_frequency_hz=best.frequency_hz,
        argmax_overlap=float(modal_assurance(reference, best.right[positions])),
    )


def band_worst(spectrum: Spectrum) -> float:
    """Worst real part over the WHOLE frozen band, with no shape criterion.

    The independent check H2C is asked to report beside the family envelope. It
    depends on no tracking decision at all.
    """

    candidates = band_candidates(spectrum, BAND_HZ)
    return float(max(m.real for m in candidates)) if candidates else float("nan")
