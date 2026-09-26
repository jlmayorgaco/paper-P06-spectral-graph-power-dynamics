"""Write hashes for preregistration inputs after the audit is frozen."""
from __future__ import annotations

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
PREREG = CAMPAIGN / "prereg"
files = [
    PREREG / "policies.json",
    PREREG / "blind_predictions.json",
    PREREG / "uncertainty_weights_manifest.json",
    PREREG / "topology_holdout.json",
    CAMPAIGN / "docs" / "PREREGISTRATION.md",
]
lines = []
for path in files:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    lines.append(f"{digest}  {path.relative_to(CAMPAIGN).as_posix()}")
(PREREG / "sha256.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"hashed {len(files)} preregistration inputs")
