"""Shared definitions of the planning/design campaign (PD).

Preregistration: docs/20260913_PLANNING_DESIGN_PREREG.md
Frozen rules:    configs/ias2026/planning_design_prereg_v1.yaml

PD1 re-runs the E34 M1 converter-only retune per operating point at the 240
unstable samples of E35; PD2 re-runs E9 (plan-level vs single-boundary design)
at census scale with a larger budget and a globalized step rule.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
ROOT = EXPERIMENTS.parent
for p in (EXPERIMENTS, EXPERIMENTS / "post_cumulant_validation", EXPERIMENTS / "connected_cumulants"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

RESULTS = ROOT / "results" / "PLANNING_DESIGN"
E35_DIR = ROOT / "outputs/ias2026/final_validation_overnight_20260910T003225/E35_MC_operating"
E35_TABLE = E35_DIR / "E35_MC_operating_2000.parquet"
CDW_EXPERIMENTS = ROOT / "contextual_dynamic_weakness" / "experiments" / "cdw"


def out_dir(name: str) -> Path:
    path = RESULTS / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def load_config() -> dict:
    import yaml

    return yaml.safe_load((ROOT / "configs/ias2026/planning_design_prereg_v1.yaml").read_text(encoding="utf-8"))
