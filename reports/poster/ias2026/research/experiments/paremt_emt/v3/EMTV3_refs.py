# ruff: noqa: E501  -- test tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - V3 stopped at gate E0 (commit 229613de). Any use needs a new preregistration.
"""V3 references for the new blind SG holdout B-SG38 (prereg V3 section 5; .venv/tx3-analysis).

Uses the V2 reference functions unchanged (experiments/paremt_emt/v2/EMTV2_refs.py: quasi-static and
dynamic-line networks, canonical classes, Radau). Bus 38, S = canonical base-PF dispatch.
Writes results/EMTV3/refs/B-SG38_{qs,dyn}.npz and refs_manifest.json.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v2"))

import EMTV2_refs as R2  # noqa: E402
import numpy as np  # noqa: E402

OUT = R2.RESEARCH / "results" / "EMTV3" / "refs"
SG_CASES = {"B-SG38": (38, 7.64783381 + 1.23275808j)}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    man = {"files": {}, "meta": {}}
    for cid, (bus, s) in SG_CASES.items():
        qs, dyn, meta = R2.sg_refs(bus, s)
        for kind, d in (("qs", qs), ("dyn", dyn)):
            p = OUT / f"{cid}_{kind}.npz"
            np.savez_compressed(p, **d)
            man["files"][str(p.relative_to(R2.RESEARCH)).replace("\\", "/")] = hashlib.sha256(p.read_bytes()).hexdigest()
        man["meta"][cid] = meta
        print(cid, "done")
    (OUT / "refs_manifest.json").write_text(json.dumps(man, indent=1, default=float), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
