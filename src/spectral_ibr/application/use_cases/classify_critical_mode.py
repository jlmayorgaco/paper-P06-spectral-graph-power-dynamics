from spectral_ibr.domain.entities.mode import Mode
from spectral_ibr.domain.entities.pole import Pole
from spectral_ibr.domain.services.mode_classifier import classify_by_control_participation


def classify_critical_mode(pole: Pole, pi_c: float) -> Mode:
    return Mode(
        pole=pole,
        family=classify_by_control_participation(pi_c),
        pi_c=pi_c,
    )

