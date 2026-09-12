# ruff: noqa: E501  -- report tables kept on one line
"""Determinism record for V2: compare a snapshot of results/EMTV2 (first execution) with the rerun.

CSV files: byte comparison. JSON files: comparison after removing provenance fields (git HEAD,
hashes, wall times). NPZ files: array equality per key. Usage: EMTV2_identity.py <snapshot_dir>.
Writes results/EMTV2/EMTV2_rerun_identity.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parents[2] / "results" / "EMTV2"
VOLATILE = {"provenance", "git_head", "paremt_working_copy_diff_sha", "environment_hash", "files"}


def strip(obj):
    if isinstance(obj, dict):
        return {k: strip(v) for k, v in obj.items() if k not in VOLATILE and not k.endswith("_wall_s")}
    if isinstance(obj, list):
        return [strip(v) for v in obj]
    return obj


def main() -> int:
    snap = Path(sys.argv[1])
    out = {}
    for p in sorted(snap.rglob("*")):
        if p.is_dir():
            continue
        rel = p.relative_to(snap)
        q = RES / rel
        key = str(rel).replace("\\", "/")
        if not q.exists():
            out[key] = "MISSING in rerun"
            continue
        if p.suffix == ".npz":
            a, b = np.load(p), np.load(q)
            out[key] = "identical arrays" if sorted(a.files) == sorted(b.files) and all(np.array_equal(a[k], b[k]) for k in a.files) else "DIFFERENT"
        elif p.suffix == ".json":
            out[key] = "identical (provenance excluded)" if strip(json.loads(p.read_text())) == strip(json.loads(q.read_text())) else "DIFFERENT"
        else:
            out[key] = "byte-identical" if p.read_bytes() == q.read_bytes() else "DIFFERENT"
    (RES / "EMTV2_rerun_identity.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
