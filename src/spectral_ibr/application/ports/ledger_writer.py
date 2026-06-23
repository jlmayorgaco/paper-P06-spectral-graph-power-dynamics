from pathlib import Path
from typing import Protocol


class ILedgerWriter(Protocol):
    def write(self, records: list[dict[str, object]], destination: Path) -> None:
        """Persist experiment records."""

