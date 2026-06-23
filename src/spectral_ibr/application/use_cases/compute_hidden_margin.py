from spectral_ibr.domain.certificates.hidden_margin import HiddenMarginCertificate


def compute_hidden_margin(
    mode_id: str,
    observed_margin: float,
    hidden_margin: float,
) -> HiddenMarginCertificate:
    return HiddenMarginCertificate(mode_id, observed_margin, hidden_margin)

