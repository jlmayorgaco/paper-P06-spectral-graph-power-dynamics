"""Capture additional transitive source inputs without changing the frozen manifest."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SNAP = HERE / "source_snapshot"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    assert not SNAP.exists(), "source snapshot already captured; refusing to overwrite"
    manifest = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    rels = set(manifest["model_sources_sha256"])
    for dirname in ("src/pd39", "src/bnd_design_e", "src/bnd_model_expN"):
        rels.update(p.relative_to(ROOT).as_posix() for p in (ROOT / dirname).rglob("*.jl"))
    rels.update(("Project.toml", "Manifest.toml"))
    report = {}
    for rel in sorted(rels):
        src = ROOT / rel
        assert src.is_file(), rel
        dst = SNAP / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src,dst)
        report[rel] = {"sha256":digest(src),"bytes":src.stat().st_size}
    (HERE / "POSTRUN_SOURCE_SNAPSHOT.json").write_text(
        json.dumps({"scope":"post-run transitive source snapshot; frozen preregistration unchanged",
                    "files":report},indent=2)+"\n",encoding="utf-8")
    print("SNAPSHOT_DONE",len(report),"files")


if __name__ == "__main__":
    main()
