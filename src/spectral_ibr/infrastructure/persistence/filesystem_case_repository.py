from pathlib import Path

from spectral_ibr.domain.entities.linearized_model import LinearizedModel
from spectral_ibr.domain.entities.power_system_case import PowerSystemCase


class FilesystemCaseRepository:
    def __init__(self, root: Path) -> None:
        self.root = root

    def load(self, case: PowerSystemCase) -> LinearizedModel:
        raise NotImplementedError(f"Load case {case.name} from {self.root}.")

