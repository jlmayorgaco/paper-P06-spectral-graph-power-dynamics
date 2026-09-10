"""Run directory for the journal gates, after tag IAS2026_TRACKA_F8_F12_POST_F7_FREEZE.

The post-F7 run directory is frozen with that tag; the gates write their
manifests under outputs/ias2026_journal_gates/ instead.
"""

from __future__ import annotations

import time
from pathlib import Path

from _overnight import RESEARCH, Experiment

GATE_ROOT = RESEARCH / "outputs" / "ias2026_journal_gates"


def gate_run_directory() -> Path:
    GATE_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = GATE_ROOT / "CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = GATE_ROOT / f"gates_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


class GateExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = gate_run_directory()
        self.directory = self.root / self.name
        self.directory.mkdir(parents=True, exist_ok=True)
