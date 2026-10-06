"""Freeze provenance and SHA-256 of every delivered closure artifact."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parent.parent


def cmd(*args: str) -> str:
    return subprocess.check_output(args, cwd=repo, text=True, encoding="utf-8").strip()


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    artifacts = {
        str(path.relative_to(here)).replace("\\", "/"): sha(path)
        for path in sorted(here.rglob("*"))
        if path.is_file() and path.name != "CLOSURE_MANIFEST.json"
        and "__pycache__" not in path.parts
    }
    status = cmd("git", "status", "--porcelain=v1")
    packages = {}
    for name in ("numpy", "scipy", "matplotlib", "pandas"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "not installed"
    manifest = {
        "experiment": str(here.relative_to(repo)).replace("\\", "/"),
        "parent_experiment": "experiments/latency_robust_pll_codesign_20261003",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": cmd("git", "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "git_status_sha256": hashlib.sha256(status.encode()).hexdigest(),
        "julia_version": cmd("julia", "--version"),
        "python_packages": packages,
        "artifact_hashes_sha256": artifacts,
    }
    (here / "CLOSURE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"manifest artifacts={len(artifacts)} head={manifest['git_head']}")


if __name__ == "__main__":
    main()
