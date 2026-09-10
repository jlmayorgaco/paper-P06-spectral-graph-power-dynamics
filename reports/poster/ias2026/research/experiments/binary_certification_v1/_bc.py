"""Shared paths and run manifests for the binary-certification programme (BC00-BC12).

Frozen F- and G-results are read, never written. Every BC script writes to
results/BC/<ID>/ and its manifest to outputs/binary_certification_v1/<run>/<ID>/.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import _bootstrap  # noqa: E402,F401  (puts src/ on the path)
from _overnight import RESEARCH, Experiment  # noqa: E402

ROOT = RESEARCH
RESULTS = ROOT / "results" / "BC"
CONFIGS = ROOT / "configs" / "binary_certification_v1"
BC_ROOT = ROOT / "outputs" / "binary_certification_v1"
WORKERS = 12  # leave cores for other jobs on the machine; never stop foreign jobs


def run_directory() -> Path:
    BC_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = BC_ROOT / "CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = BC_ROOT / f"bc_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


class BCExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = run_directory()
        self.directory = self.root / self.name
        self.directory.mkdir(parents=True, exist_ok=True)


def out_dir(bc_id: str) -> Path:
    path = RESULTS / bc_id
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def load_yaml(path: Path) -> dict:
    import yaml

    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))
