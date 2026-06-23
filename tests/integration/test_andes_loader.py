from spectral_ibr.infrastructure.andes.andes_case_loader import AndesCaseLoader
from spectral_ibr.infrastructure.andes.sg_to_ibr_replacer import SgToIbrReplacer


def test_andes_loader_adapter_exists() -> None:
    assert AndesCaseLoader() is not None


def test_sg_to_ibr_replacer_adapter_exists() -> None:
    assert SgToIbrReplacer() is not None
