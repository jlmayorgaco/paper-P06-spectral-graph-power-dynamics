"""Record analytical LP prediction versus exact spectral correction for every trial."""
from __future__ import annotations

import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def first(path: Path):
    with path.open(newline="") as stream:
        return next(csv.DictReader(stream))


out = []
for prediction in sorted(HERE.glob("PREDICTION_*.csv")):
    p = first(prediction)
    cid = p["candidate_id"]
    spec = HERE / "evaluations" / cid / "SPECTRAL.csv"
    if not spec.exists():
        continue
    roots = HERE / "evaluations" / cid / "ROOTS.csv"
    with roots.open(newline="") as stream:
        rr = list(csv.DictReader(stream))
    corrected = min(float(r["local_crossing_ms"]) for r in rr)
    predicted = float(p["predicted_tau_ms"])
    validation = HERE / "event_validations" / cid / "Q0_RESULT.toml"
    out.append(dict(candidate_id=cid, parent=p["parent"], predicted_tau_ms=predicted,
                    first_refined_crossing_ms=corrected,
                    prediction_error_ms=corrected-predicted,
                    root_count_in_catalogue=len(rr),
                    five_event_record_present=validation.exists(),
                    status="NUMERICAL_CORRECTION_NOT_GLOBAL_CERTIFICATE"))
with (HERE / "TABLE_F12_PREDICTOR_CORRECTION.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
print("correction audit", len(out), "trials")
