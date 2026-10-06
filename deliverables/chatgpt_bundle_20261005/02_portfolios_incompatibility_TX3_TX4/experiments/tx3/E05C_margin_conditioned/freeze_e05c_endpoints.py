from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05C_margin_conditioned"
PREREG = ARTIFACT / "preregistration"
LIMITS = ARTIFACT / "baseline_calibration" / "E05C_BASELINE_STRESS_LIMITS.parquet"
TRACE = ARTIFACT / "baseline_calibration" / "E05C_BASELINE_CALIBRATION_TRACE.parquet"


def sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    if (ARTIFACT / "holdout" / "e05c_connected_pole_externalities.parquet").exists():
        raise RuntimeError("cannot freeze stress endpoints after coalition results exist")
    limits = pd.read_parquet(LIMITS).sort_values("seed_id")
    if len(limits) != 8 or not limits["status"].eq("SUCCESS").all():
        raise RuntimeError("eight successful baseline-only limits are required")
    records = limits.to_dict("records")
    payload = {
        "freeze_id": "TX3-E05C-STRESS-ENDPOINTS-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "calibration_population": "EMPTY coalition only",
        "coalition_outcomes_inspected_before_freeze": False,
        "limits_path": str(LIMITS.relative_to(ROOT)).replace("\\", "/"),
        "limits_sha256": sha256(LIMITS),
        "trace_path": str(TRACE.relative_to(ROOT)).replace("\\", "/"),
        "trace_sha256": sha256(TRACE),
        "seed_count": 8,
        "gamma_end_range": [float(limits["gamma_end"].min()), float(limits["gamma_end"].max())],
        "endpoint_reason_counts": limits["endpoint_reason"].value_counts().to_dict(),
        "limits": records,
    }
    output = PREREG / "E05C_STRESS_ENDPOINT_FREEZE.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
