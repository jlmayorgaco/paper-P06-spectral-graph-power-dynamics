from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

_TRACKED = ("numpy", "scipy", "pandas", "matplotlib", "networkx", "pyarrow")


def _git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _versions() -> dict[str, str]:
    found: dict[str, str] = {}
    for name in _TRACKED:
        try:
            found[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            found[name] = "absent"
    return found


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass
class Manifest:
    """Everything needed to re-run one experiment and get the same numbers."""

    experiment: str
    seed: int
    config: dict[str, Any]
    python: str = field(default_factory=lambda: sys.version.split()[0])
    platform_name: str = field(default_factory=platform.platform)
    git_commit: str = field(default_factory=_git_commit)
    dependencies: dict[str, str] = field(default_factory=_versions)
    started_utc: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    results: dict[str, Any] = field(default_factory=dict)
    status: str = "RUNNING"

    @property
    def config_sha256(self) -> str:
        return config_hash(self.config)

    def finish(self, status: str, **results: Any) -> Manifest:
        self.status = status
        self.results.update(results)
        return self

    def write(self, directory: Path) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        payload["config_sha256"] = self.config_sha256
        payload["finished_utc"] = datetime.now(UTC).isoformat()
        path = directory / f"{self.experiment}.json"
        path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
        return path
