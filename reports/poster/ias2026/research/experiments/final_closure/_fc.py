"""Shared paths for the final closure campaign (FC00-FC20).

Everything new is written under outputs/ias2026/final_math_nonlinear_validation_<ts>/.
Named deliverables requested by the campaign brief are also copied to
results/ (append-only). Frozen F/G/E results are read, never written.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
for p in (EXPERIMENTS, EXPERIMENTS / "binary_certification_v1"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import _bootstrap  # noqa: E402,F401  (puts src/ on the path)
from _overnight import RESEARCH, Experiment  # noqa: E402

ROOT = RESEARCH
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs" / "ias2026"
OUT_ROOT = ROOT / "outputs" / "ias2026"
WORKERS = int(__import__("os").environ.get("FC_WORKERS", "16"))
CORE = (30, 33, 35, 37)


def run_directory() -> Path:
    stamp = OUT_ROOT / "FINAL_CLOSURE_CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = OUT_ROOT / f"final_math_nonlinear_validation_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


RUN = run_directory()


class FCExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = RUN
        self.directory = RUN / self.name
        self.directory.mkdir(parents=True, exist_ok=True)


def out_dir(name: str) -> Path:
    path = RUN / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def publish(path: Path, name: str | None = None) -> Path:
    """Copy a campaign output to results/ under its declared deliverable name."""

    target = RESULTS / (name or Path(path).name)
    shutil.copy2(path, target)
    return target


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def label(members) -> str:
    return "+".join(map(str, sorted(members))) or "BASE"
