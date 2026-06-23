import csv
from pathlib import Path


class CsvLedgerWriter:
    def write(self, records: list[dict[str, object]], destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not records:
            destination.write_text("", encoding="utf-8")
            return
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

