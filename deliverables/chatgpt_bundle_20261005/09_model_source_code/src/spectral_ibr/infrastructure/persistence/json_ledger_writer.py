import json
from pathlib import Path


class JsonLedgerWriter:
    def write(self, records: list[dict[str, object]], destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(records, indent=2), encoding="utf-8")

