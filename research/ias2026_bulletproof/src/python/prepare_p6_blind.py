"""Freeze P6 discovery/holdout predictions before revealing holdout labels."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


CAMPAIGN = Path(__file__).resolve().parents[2]
P2 = CAMPAIGN / "raw" / "p2" / "p2_powerdynamics_portfolios.csv"
PREREG = CAMPAIGN / "prereg"
RAW = CAMPAIGN / "raw" / "p6"
PREREG.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rows = sorted(csv.DictReader(P2.open(encoding="utf-8", newline="")), key=lambda r: r["portfolio"])
    discovery = rows[:12]
    holdout = rows[12:]
    # The only labels used to construct the frozen predictor are the discovery
    # labels. Holdout labels are deliberately not serialized into this file.
    train_y = np.array([r["stable"].lower() == "true" for r in discovery], dtype=float)
    majority_stable = bool(train_y.mean() >= 0.5)
    x = np.array([int(r["replaced_count"]) for r in discovery], dtype=float)
    alpha = np.array([float(r["alpha_transverse"]) for r in discovery], dtype=float)
    freq = np.array([float(r["critical_frequency_hz"]) for r in discovery], dtype=float)
    alpha_coeff = np.polyfit(x, alpha, 1).tolist() if len(set(x)) > 1 else [0.0, float(alpha.mean())]
    freq_coeff = np.polyfit(x, freq, 1).tolist() if len(set(x)) > 1 else [0.0, float(freq.mean())]
    predictions = {}
    for row in rows:
        count = int(row["replaced_count"])
        predictions[row["portfolio"]] = {
            "split": "discovery" if row in discovery else "holdout",
            "predicted_stable": majority_stable,
            "predicted_alpha_transverse": float(np.polyval(alpha_coeff, count)),
            "predicted_critical_frequency_hz": float(np.polyval(freq_coeff, count)),
        }
    payload = {
        "protocol": "P6 blind holdout v1",
        "seed": 20260919,
        "ordering": "lexicographic portfolio key",
        "discovery_portfolios": [r["portfolio"] for r in discovery],
        "holdout_portfolios": [r["portfolio"] for r in holdout],
        "prediction_rule": "majority stable label on discovery split; linear count-only roots fitted on discovery split",
        "known_branch_repairs_excluded": [0, 1, 13, 43],
        "predictions": predictions,
    }
    path = PREREG / "blind_predictions_P6.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (RAW / "p6_blind_predictions.sha256").write_text(f"{digest}  prereg/blind_predictions_P6.json\n", encoding="utf-8")
    (RAW / "p6_blind_manifest.json").write_text(
        json.dumps({
            "prediction_path": "prereg/blind_predictions_P6.json",
            "sha256": digest,
            "discovery_count": len(discovery),
            "holdout_count": len(holdout),
            "status": "FROZEN_BEFORE_HOLDOUT_REVEAL",
            "git_commit": "PENDING_COMMIT",
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"P6_PREDICTIONS_FROZEN sha256={digest} discovery={len(discovery)} holdout={len(holdout)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
