"""Registered H0 blocker-antichain definitions for the TX4 audit."""

from __future__ import annotations

import math
from typing import Mapping


H4_LABEL = "30+33+35+37"


def portfolio_members(label: str) -> frozenset[int]:
    """Return the bus set encoded by a master-table portfolio label."""
    if label in ("", "BASE"):
        return frozenset()
    return frozenset(int(part) for part in label.split("+"))


def _sort_key(label: str) -> tuple[int, str]:
    return (len(portfolio_members(label)), label)


def classify_blockers(
    alpha_em: Mapping[str, float | int | None],
    *,
    threshold: float = 0.0,
    h4_label: str = H4_LABEL,
) -> dict[str, object]:
    """Classify U0 and its inclusion-minimal blocker antichain H0.

    Missing, non-finite, or non-numeric alpha values are not members of U0;
    the caller retains the associated row status separately. The registered
    blocker convention uses ``alpha_EM >= 0`` exactly, so the default
    threshold is zero rather than the legacy 1e-8 engineering flag threshold.
    """
    u0: list[str] = []
    for label, value in alpha_em.items():
        try:
            numeric = float(value) if value is not None else float("nan")
        except (TypeError, ValueError):
            numeric = float("nan")
        if math.isfinite(numeric) and numeric >= threshold:
            u0.append(label)

    u0 = sorted(set(u0), key=_sort_key)
    u0_sets = {label: portfolio_members(label) for label in u0}
    h0 = [
        label
        for label in u0
        if not any(
            other != label
            and other_set < u0_sets[label]
            for other, other_set in u0_sets.items()
        )
    ]
    h0 = sorted(h0, key=_sort_key)
    h4_unstable = h4_label in u0
    h4_present = h4_label in h0
    exact_h4 = h0 == [h4_label]
    noncomposable = any(len(portfolio_members(label)) >= 2 for label in h0)
    kappa = min((len(portfolio_members(label)) for label in h0), default=None)
    result = {
        "U0": u0,
        "H0": h0,
        "H4_UNSTABLE": h4_unstable,
        "H4_PRESENT": h4_present,
        "EXACT_H4": exact_h4,
        "NONCOMPOSABLE": noncomposable,
        "KAPPA": kappa,
    }
    if exact_h4:
        assert h4_present and h4_unstable and noncomposable
    return result

