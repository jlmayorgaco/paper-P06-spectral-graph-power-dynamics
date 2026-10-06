"""Capture the operating-point and network CSV inputs used by the model."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SNAP = HERE / "input_snapshot"
INPUTS = [
    "reports/experiment_D/inputs/bus.csv",
    "reports/experiment_D/inputs/branch.csv",
    "reports/experiment_D/inputs/load.csv",
    "reports/experiment_D/inputs/machine.csv",
    "reports/experiment_A/matrices/bus33_baseline_equilibrium.csv",
    "reports/experiment_N/TABLE_N01_original_operating_point.csv",
]


def main() -> None:
    assert not SNAP.exists(), "input snapshot already exists"
    out = {}
    for rel in INPUTS:
        src = ROOT / rel
        assert src.is_file(), rel
        dest = SNAP / rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest)
        out[rel] = {"sha256":hashlib.sha256(src.read_bytes()).hexdigest(),
                    "bytes":src.stat().st_size}
    (HERE/"POSTRUN_INPUT_SNAPSHOT.json").write_text(json.dumps({
        "scope":"post-run model CSV inputs; not a rewrite of the frozen preregistration",
        "files":out},indent=2)+"\n",encoding="utf-8")
    print("INPUT_SNAPSHOT_DONE",len(out),"files")


if __name__ == "__main__":
    main()
