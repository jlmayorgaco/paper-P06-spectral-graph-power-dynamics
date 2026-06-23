from spectral_ibr.domain.certificates.braess_resonant import BraessResonantCertificate


def audit_line_reinforcement(
    action_id: str,
    sensitivity: float,
    resonance_gap: float,
    threshold: float,
) -> BraessResonantCertificate:
    return BraessResonantCertificate(action_id, sensitivity, resonance_gap, threshold)

