"""Shared paths for the TX4 paper Monte Carlo validation (MC01-MC06).

Everything is written under outputs/ias2026/paper_mc_<ts>/ (stamp file
PAPER_MC_CURRENT_RUN). The frozen pass rules are
configs/ias2026/paper_mc_validation_v1.yaml; results never change them.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
for p in (
    EXPERIMENTS,
    EXPERIMENTS / "binary_certification_v1",
    EXPERIMENTS / "final_closure",
):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import _bootstrap  # noqa: E402,F401  (puts src/ on the path)
from _overnight import RESEARCH, Experiment  # noqa: E402

ROOT = RESEARCH
CONFIGS = ROOT / "configs" / "ias2026"
OUT_ROOT = ROOT / "outputs" / "ias2026"
WORKERS = int(os.environ.get("MC_WORKERS", "16"))
CFG = yaml.safe_load(
    (CONFIGS / "paper_mc_validation_v1.yaml").read_text(encoding="utf-8")
)
SEED = int(CFG["seed"])
OMEGA_B = float(CFG["omega_b"])


def run_directory() -> Path:
    stamp = OUT_ROOT / "PAPER_MC_CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = OUT_ROOT / f"paper_mc_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


RUN = run_directory()


class MCExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = RUN
        self.directory = RUN / self.name
        self.directory.mkdir(parents=True, exist_ok=True)


def out_dir(name: str) -> Path:
    path = RUN / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")
