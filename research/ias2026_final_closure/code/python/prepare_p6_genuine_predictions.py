"""Freeze genuine P6 holdout predictions after discovery labels and before reveal."""

from __future__ import annotations

import csv
import hashlib
import json

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
inputs = CAMPAIGN / "prereg" / "p6_genuine_holdout_inputs.csv"
discovery = CAMPAIGN / "raw" / "p6" / "p6_genuine_discovery_labels.csv"
out = CAMPAIGN / "prereg" / "p6_genuine_predictions.json"

input_digest = hashlib.sha256(inputs.read_bytes()).hexdigest()
rows = list(csv.DictReader(discovery.open(encoding="utf-8", newline="")))
majority_stable = sum(r["stable"].lower() == "true" for r in rows) >= (len(rows) / 2)
holdout = [r for r in csv.DictReader(inputs.open(encoding="utf-8", newline="")) if r["phase"] == "holdout"]
payload = {
    "protocol": "P6 genuine holdout v2",
    "model_scope": "official PowerDynamics GFL11 alternative-model harness",
    "input_sha256": input_digest,
    "discovery_case_ids": [r["case_id"] for r in rows],
    "holdout_case_ids": [r["case_id"] for r in holdout],
    "prediction_rule": "majority stable label on discovery phase only",
    "predictions": {r["case_id"]: {"predicted_stable": majority_stable, "portfolio": r["portfolio"], "gfl_gain": float(r["gfl_gain"])} for r in holdout},
}
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
digest = hashlib.sha256(out.read_bytes()).hexdigest()
(CAMPAIGN / "raw" / "p6" / "p6_genuine_prediction_sha256.txt").write_text(f"{digest}  prereg/p6_genuine_predictions.json\n", encoding="utf-8")
print(f"P6_GENUINE_PREDICTIONS_FROZEN sha256={digest} holdout={len(holdout)}")
