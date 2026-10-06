"""Hash all experiment artifacts, excluding this generated manifest and ZIPs."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
records = []
for path in sorted(HERE.rglob("*")):
    if not path.is_file() or path.name == "RESULT_HASHES.json" or path.suffix.lower() == ".zip":
        continue
    data = path.read_bytes()
    records.append({
        "path": path.relative_to(HERE).as_posix(),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    })
out = {"scope": "new Mega experiment directory only", "file_count": len(records), "files": records}
(HERE / "RESULT_HASHES.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
print("HASHED", len(records), "files")
