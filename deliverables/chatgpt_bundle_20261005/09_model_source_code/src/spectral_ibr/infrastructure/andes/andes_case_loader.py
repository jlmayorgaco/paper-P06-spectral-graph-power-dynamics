from pathlib import Path

from spectral_ibr.domain.entities.linearized_model import LinearizedModel


class AndesCaseLoader:
    def load_xlsx(self, path: Path) -> LinearizedModel:
        raise NotImplementedError(f"Wire ANDES case loading for {path}.")

