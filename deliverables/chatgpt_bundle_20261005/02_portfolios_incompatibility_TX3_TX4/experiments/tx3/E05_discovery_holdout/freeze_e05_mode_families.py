from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05_mechanism_consequence"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    holdout_outputs = list((ARTIFACT / "holdout").glob("*.parquet"))
    if holdout_outputs:
        raise RuntimeError("mode-family freeze cannot be regenerated after holdout execution")
    numerical = json.loads((ARTIFACT / "preregistration" / "E05_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    candidates = json.loads((ARTIFACT / "preregistration" / "E05_CANDIDATE_FREEZE.json").read_text(encoding="utf-8"))
    catalog = json.loads((ARTIFACT / "discovery" / "REFERENCE_MODE_CATALOG.json").read_text(encoding="utf-8"))
    consequence_path = ARTIFACT / "discovery" / "mode_consequences_all_families.parquet"
    frame = pd.read_parquet(consequence_path)
    selected = []
    minimum_bmac = float(numerical["minimum_tracking_bmac"])
    minimum_coverage = int(numerical["discovery_mode_family_minimum_valid_points"])
    catalog_by_id = {record["reference_mode_id"]: record for record in catalog["modes"]}
    for candidate in candidates["candidates"]:
        candidate_frame = frame[frame["candidate_id"] == candidate["candidate_id"]]
        ranking = []
        for mode_id, group in candidate_frame.groupby("reference_mode_id"):
            valid = group[group["minimum_tracking_bmac"] >= minimum_bmac]
            if len(valid) < minimum_coverage:
                continue
            ranking.append(
                {
                    "reference_mode_id": mode_id,
                    "valid_discovery_points": len(valid),
                    "median_absolute_damping_externality": float(valid["damping_externality"].abs().median()),
                    "median_damping_externality": float(valid["damping_externality"].median()),
                    "median_real_part_externality_per_s": float(valid["real_part_externality_per_s"].median()),
                    "minimum_tracking_bmac": float(valid["minimum_tracking_bmac"].min()),
                }
            )
        if not ranking:
            raise RuntimeError(f"no discovery mode family satisfies tracking for {candidate['candidate_id']}")
        winner = max(ranking, key=lambda record: (record["median_absolute_damping_externality"], record["reference_mode_id"]))
        winner["candidate_id"] = candidate["candidate_id"]
        winner["coalition_key"] = candidate["coalition_key"]
        winner["selection_rule"] = "largest discovery median absolute damping-ratio Mobius externality among reference mode families with >=12/16 valid BMAC-tracked points"
        winner["reference_mode"] = catalog_by_id[winner["reference_mode_id"]]
        selected.append(winner)
    payload = {
        "freeze_id": "TX3-E05-MODE-FAMILIES-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "holdout_consequences_inspected_before_freeze": False,
        "discovery_consequence_sha256": sha256(consequence_path),
        "minimum_tracking_bmac": minimum_bmac,
        "selected_mode_families": selected,
    }
    path = ARTIFACT / "preregistration" / "E05_MODE_FAMILY_FREEZE.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
