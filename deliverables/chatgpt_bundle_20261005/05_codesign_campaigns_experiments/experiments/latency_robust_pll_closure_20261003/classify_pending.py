"""Gate the stored 41 ms predictor using complete roots and five events."""
from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "latency_robust_pll_codesign_20261003"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    cid = "pending_41ms"
    evaldir = HERE / "evaluations" / cid
    spectral = rows(evaldir / "SPECTRAL.csv")[0]
    coverage = rows(evaldir / "ROOT_COVERAGE.csv")[0]
    modes = rows(evaldir / "ROOTS.csv")
    result = tomllib.loads((HERE / "event_validations" / cid / "Q0_RESULT.toml").read_text())
    old = json.loads((OLD / "RESULT_SUMMARY.json").read_text())
    manifest = json.loads((HERE / "PARENT_EXPERIMENT_MANIFEST.json").read_text())
    tau = min(float(m["local_crossing_ms"]) for m in modes)
    events = result["events"]
    gain_bounds = bool(result["gain_bounds_pass"])
    spectrum = (coverage["status"] == "COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
                and spectral["complete_contour_bracket_pass"].lower() == "true"
                and float(spectral["alpha_zero_delay_s_inv"]) <= -0.05)
    actuator = min(float(e["min_SG_actuator_fraction_slack"]) for e in events)
    physical = all(e["pass"] for e in events)
    if not gain_bounds:
        status = "REJECTED_GAIN_BOUND"
    elif not spectrum:
        status = "REJECTED_SPECTRAL"
    elif not physical and actuator < manifest["event_limits"]["min_SG_actuator_fraction_slack"]:
        status = "REJECTED_ACTUATOR"
    elif not physical:
        status = "REJECTED_EVENT"
    elif tau <= float(old["best_local_numerical_uniform_latency_threshold_ms"]):
        status = "REJECTED_SPECTRAL"
    else:
        status = "ACCEPTED"
    record = dict(candidate_id=cid, status=status, true_tau_crit_ms=tau,
                  previous_best_tau_ms=old["best_local_numerical_uniform_latency_threshold_ms"],
                  improvement_ms=tau-float(old["best_local_numerical_uniform_latency_threshold_ms"]),
                  alpha_zero_delay_s_inv=spectral["alpha_zero_delay_s_inv"],
                  complete_root_coverage=coverage["status"],
                  modal_pair_count=len(modes),
                  first_four_crossing_ms=";".join(m["local_crossing_ms"] for m in modes[:4]),
                  all_five_events_pass=physical,
                  actuator_slack=actuator,
                  actuator_limit=manifest["event_limits"]["min_SG_actuator_fraction_slack"],
                  max_frequency_deviation_hz=max(float(e["F_peak_Hz"]) for e in events),
                  max_rocof_hz_s=max(float(e["RoCoF_peak_Hz_s"]) for e in events),
                  gfl_percent=result["GFL_percent"],gfl_mw=result["GFL_MW"],
                  retained_sg_mw=result["retained_SG_MW"])
    for name in ("F01_PENDING_CANDIDATE_VALIDATION.csv", "TABLE_F1_PENDING_CANDIDATE.csv"):
        with (HERE / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(record))
            writer.writeheader()
            writer.writerow(record)
    print(record)


if __name__ == "__main__":
    main()
