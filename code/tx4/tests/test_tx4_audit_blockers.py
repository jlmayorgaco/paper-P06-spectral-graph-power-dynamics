from __future__ import annotations

import math

from tx4_audit_blockers import classify_blockers


H4 = "30+33+35+37"


def test_empty_u0_has_no_minimal_blocker():
    result = classify_blockers({"BASE": -0.1, "30": -0.2, H4: -0.3})
    assert result["U0"] == []
    assert result["H0"] == []
    assert result["H4_UNSTABLE"] is False
    assert result["H4_PRESENT"] is False
    assert result["EXACT_H4"] is False
    assert result["NONCOMPOSABLE"] is False
    assert result["KAPPA"] is None


def test_singleton_h4_is_exact_h4_and_noncomposable_by_cardinality():
    result = classify_blockers({"BASE": -0.1, "30": -0.2, H4: 0.0})
    assert result["U0"] == [H4]
    assert result["H0"] == [H4]
    assert result["H4_UNSTABLE"] is True
    assert result["H4_PRESENT"] is True
    assert result["EXACT_H4"] is True
    assert result["NONCOMPOSABLE"] is True
    assert result["KAPPA"] == 4


def test_proper_singleton_is_not_h4_present():
    result = classify_blockers({"BASE": -0.1, "30": 0.0, H4: -0.2})
    assert result["H0"] == ["30"]
    assert result["H4_UNSTABLE"] is False
    assert result["H4_PRESENT"] is False
    assert result["EXACT_H4"] is False
    assert result["NONCOMPOSABLE"] is False
    assert result["KAPPA"] == 1


def test_composite_blocker_is_noncomposable():
    result = classify_blockers({"BASE": -0.1, "30+33": 0.2, H4: 0.3})
    assert result["H0"] == ["30+33"]
    assert result["H4_UNSTABLE"] is True
    assert result["H4_PRESENT"] is False
    assert result["EXACT_H4"] is False
    assert result["NONCOMPOSABLE"] is True
    assert result["KAPPA"] == 2


def test_mixed_minimal_antichain_preserves_both_blockers():
    result = classify_blockers({"30": 0.1, "33+35": 0.2, H4: 0.3})
    assert result["H0"] == ["30", "33+35"]
    assert result["H4_UNSTABLE"] is True
    assert result["H4_PRESENT"] is False
    assert result["EXACT_H4"] is False
    assert result["NONCOMPOSABLE"] is True
    assert result["KAPPA"] == 1


def test_incomparable_composite_blockers_are_both_minimal():
    result = classify_blockers({"30+33": 0.1, "30+35+37": 0.2, H4: -0.1})
    assert result["H0"] == ["30+33", "30+35+37"]
    assert result["H4_UNSTABLE"] is False
    assert result["H4_PRESENT"] is False
    assert result["EXACT_H4"] is False
    assert result["NONCOMPOSABLE"] is True
    assert result["KAPPA"] == 2
