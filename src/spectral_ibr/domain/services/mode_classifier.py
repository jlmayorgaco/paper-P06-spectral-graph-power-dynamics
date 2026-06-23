from spectral_ibr.domain.value_objects.mode_family import ModeFamily


def classify_by_control_participation(
    pi_c: float,
    threshold: float = 0.5,
) -> ModeFamily:
    return ModeFamily.CONTROL if pi_c >= threshold else ModeFamily.NETWORK

