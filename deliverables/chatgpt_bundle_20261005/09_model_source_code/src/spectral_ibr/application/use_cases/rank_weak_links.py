from spectral_ibr.domain.certificates.weak_link import WeakLinkCertificate


def rank_weak_links(certificates: list[WeakLinkCertificate]) -> list[WeakLinkCertificate]:
    return sorted(certificates, key=lambda item: item.score, reverse=True)

