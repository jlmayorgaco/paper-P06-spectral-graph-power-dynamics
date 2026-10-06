from spectral_ibr.domain.certificates.weak_node import WeakNodeCertificate


def rank_weak_nodes(certificates: list[WeakNodeCertificate]) -> list[WeakNodeCertificate]:
    return sorted(certificates, key=lambda item: item.score, reverse=True)

