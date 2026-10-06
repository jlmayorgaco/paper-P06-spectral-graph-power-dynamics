"""Hash every experiment artifact for handoff and exact local reproduction."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "EXPERIMENT_ARTIFACT_HASHES.json"
items = []
for path in sorted(p for p in HERE.rglob("*")
                   if p.is_file() and p != MANIFEST and p.suffix.lower() != ".zip"):
    rel = path.relative_to(HERE).as_posix()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    items.append({"path": rel, "size_bytes": path.stat().st_size, "sha256": digest})
payload = {
    "campaign": "delay_dressed_replacement_frontier_20261002",
    "hash_algorithm": "SHA-256",
    "hashes_include": "all files under this experiment directory except this manifest and ZIP distribution archives",
    "artifact_count": len(items),
    "artifacts": items,
}
MANIFEST.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"hashed {len(items)} experiment artifacts -> {MANIFEST.name}")
