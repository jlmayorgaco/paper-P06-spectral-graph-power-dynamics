from spectral_ibr.domain.services.mode_classifier import classify_by_control_participation
from spectral_ibr.domain.value_objects.mode_family import ModeFamily


def test_classifies_control_family_at_threshold() -> None:
    assert classify_by_control_participation(0.5) == ModeFamily.CONTROL


def test_classifies_network_family_below_threshold() -> None:
    assert classify_by_control_participation(0.49) == ModeFamily.NETWORK

