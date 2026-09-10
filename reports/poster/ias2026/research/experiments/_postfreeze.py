"""Run directory for work after the F7 freeze.

The frozen overnight run directory (and its CURRENT_RUN pointer) belong to tag
IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE and must not change. Every later
phase writes its manifests under outputs/ias2026_post_f7/ instead.
"""

from __future__ import annotations

import time
from pathlib import Path

from _overnight import RESEARCH, Experiment

POST_ROOT = RESEARCH / "outputs" / "ias2026_post_f7"


def post_run_directory() -> Path:
    POST_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = POST_ROOT / "CURRENT_RUN"
    if stamp.exists():
        path = Path(stamp.read_text(encoding="utf-8").strip())
        if path.exists():
            return path
    path = POST_ROOT / f"post_f7_{time.strftime('%Y%m%dT%H%M%S')}"
    path.mkdir(parents=True, exist_ok=True)
    stamp.write_text(str(path), encoding="utf-8")
    return path


class PostFreezeExperiment(Experiment):
    def __post_init__(self) -> None:
        self.root = post_run_directory()
        self.directory = self.root / self.name
        self.directory.mkdir(parents=True, exist_ok=True)
