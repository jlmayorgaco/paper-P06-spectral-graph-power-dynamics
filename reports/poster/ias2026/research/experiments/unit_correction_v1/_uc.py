"""Shared paths for the Phase-I unit correction (UC01-UC03).

Append-only: frozen E/F/G tables are READ, never written. Outputs go to
results/UC/<ID>/ and manifests to outputs/unit_correction_v1/<run>/<ID>/.
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
RESULTS = ROOT / "results" / "UC"
FROZEN_TABLES = ROOT / "results" / "tables"
OVERNIGHT = ROOT / "outputs" / "ias2026" / "final_validation_overnight_20260910T003225"
UC_ROOT = ROOT / "outputs" / "unit_correction_v1"
WORKERS = 12
CORE = (30, 33, 35, 37)


def run_directory() -> Path:
    UC_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = UC_ROOT / "CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = UC_ROOT / f"uc_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


class UCExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = run_directory()
        self.directory = self.root / self.name
        self.directory.mkdir(parents=True, exist_ok=True)


def out_dir(uc_id: str) -> Path:
    path = RESULTS / uc_id
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def parse(members: str) -> tuple[int, ...]:
    return (
        ()
        if members in ("BASE", "", "EMPTY")
        else tuple(int(b) for b in str(members).split("+"))
    )
