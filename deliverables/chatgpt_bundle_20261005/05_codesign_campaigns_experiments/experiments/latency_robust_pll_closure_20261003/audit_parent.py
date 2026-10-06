"""Independently reconcile the frozen predecessor's stored numerical outputs."""
from __future__ import annotations

import csv
import json
import math
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "latency_robust_pll_codesign_20261003"
SUMMARY = json.loads((OLD / "RESULT_SUMMARY.json").read_text())


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def event_metrics(design: str) -> dict:
    result = tomllib.loads((OLD / "event_validations" / design / "Q0_RESULT.toml").read_text())
    events = result["events"]
    assert len(events) == 5
    return {
        "all_five_pass": all(e["pass"] for e in events),
        "minimum_actuator_slack": min(e["min_SG_actuator_fraction_slack"] for e in events),
        "max_frequency_deviation_hz": max(e["F_peak_Hz"] for e in events),
        "max_rocof_hz_s": max(e["RoCoF_peak_Hz_s"] for e in events),
        "gfl_mw": result["GFL_MW"],
        "gfl_percent": result["GFL_percent"],
        "retained_sg_mw": result["retained_SG_MW"],
    }


def main() -> None:
    identifiers = {
        "Z": "Z_zero_delay_tuned",
        "N": "N_nominal",
        "previous_best": SUMMARY["best_fully_validated_found_id"],
    }
    records = []
    for label, design in identifiers.items():
        if label in ("Z", "N"):
            spectral = next(r for r in rows(OLD / "L0_REPRODUCTION.csv") if r["design_id"] == label)
            critical = float(spectral["tau_crit_local_ms"])
            alpha = float(spectral["alpha_zero_delay_s_inv"])
        else:
            spectral = rows(OLD / "evaluations" / design / "SPECTRAL.csv")[0]
            roots = rows(OLD / "evaluations" / design / "ROOTS.csv")
            critical = min(float(r["local_crossing_ms"]) for r in roots)
            assert math.isclose(critical, float(spectral["local_tau_crit_ms"]), abs_tol=1e-7)
            alpha = float(spectral["alpha_zero_delay_s_inv"])
        event = event_metrics(design)
        assert event["all_five_pass"]
        records.append({
            "label": label,
            "design_id": design,
            "tau_crit_ms": critical,
            "alpha_zero_delay_s_inv": alpha,
            **event,
            "status": "REPRODUCED_FROM_FROZEN_ROOTS_AND_FIVE_EVENT_RECORDS",
        })
    best = records[-1]
    for source, observed in [
        ("fixed_GFL_percent", best["gfl_percent"]),
        ("fixed_GFL_MW", best["gfl_mw"]),
        ("fixed_retained_SG_MW", best["retained_sg_mw"]),
        ("best_local_numerical_uniform_latency_threshold_ms", best["tau_crit_ms"]),
        ("best_min_SG_actuator_slack", best["minimum_actuator_slack"]),
    ]:
        assert math.isclose(float(SUMMARY[source]), observed, rel_tol=1e-8, abs_tol=1e-8), source
    with (HERE / "F00_PARENT_REPRODUCTION.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
