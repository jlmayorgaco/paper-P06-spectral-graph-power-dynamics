"""Check the internal consistency of delivered claims against generated data."""
from __future__ import annotations

import csv
import math
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    trace = rows(HERE / "TABLE_F2_FINAL_OPTIMIZATION_TRACE.csv")
    validated = [r for r in trace if r["status"] == "VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"]
    assert validated
    best = max(validated, key=lambda r: float(r["tau_crit_ms"]))
    cid = best["design_id"]
    assert (HERE / "BEST_VALIDATED_DESIGN_ID.txt").read_text().strip() == cid
    model = tomllib.loads((HERE / "designs" / f"{cid}.toml").read_text())
    event = tomllib.loads((HERE / "event_validations" / cid / "Q0_RESULT.toml").read_text())
    assert event["all_five_events_pass"] and len(event["events"]) == 5
    assert float(best["actuator_slack"]) >= 0.002
    assert all(0.0 < float(r) < 1.0 for r in model["rho"])
    for r in validated:
        d = tomllib.loads((HERE / "designs" / f"{r['design_id']}.toml").read_text())
        assert d["rho"] == model["rho"], r["design_id"]
        ev = tomllib.loads((HERE / "event_validations" / r["design_id"] / "Q0_RESULT.toml").read_text())
        assert ev["all_five_events_pass"]
        rootdir = HERE / "evaluations" / r["design_id"]
        roots = rows(rootdir / "ROOTS.csv")
        tau = min(float(x["local_crossing_ms"]) for x in roots)
        assert abs(tau - float(r["tau_crit_ms"])) < 1e-8
        cov = rootdir / "ROOT_COVERAGE.csv"
        if not cov.exists():
            cov = rootdir / "ROOT_COVERAGE_HIGH.csv"
        assert rows(cov)[0]["status"] == "COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
        assert r["boundary_contour_pass"] == "True"
    samples = [r for r in rows(HERE / "TABLE_F8_POSITIVE_DELAY_VALIDATION.csv")
               if r["design_id"] == cid and r["steps_per_delay"] == "160"]
    assert len(samples) == 3
    signs = {float(r["factor"]): r["observed_behavior"] for r in samples}
    assert signs == {0.9: "DECAY", 0.98: "DECAY", 1.02: "GROWTH"}
    rank = rows(HERE / "TABLE_F11_GAIN_ACTION_RANK.csv")
    assert rank and all(int(r["action_rank"]) <= 10 for r in rank)
    assert max(float(r["relative_action_reconstruction_error"]) for r in rank) < 1e-9
    for name in ("TABLE_F5_MULTISTART.csv", "TABLE_F6_GENERIC_OPTIMIZER_COMPARISON.csv", "TABLE_F7_PARETO_FRONTIER.csv"):
        assert rows(HERE / name)[0]["status"] == "NOT_EXECUTED"
    for name in ("FIG_F1_LATENCY_MARGIN_OPTIMIZATION.png", "FIG_F2_MODAL_FAMILY_EQUALIZATION.png",
                 "FIG_F3_NOMINAL_LATENCY_PARETO.png", "FIG_F4_ACTUATOR_LATENCY_FRONTIER.png",
                 "FIG_F5_DELAYED_TIME_DOMAIN_VALIDATION.png", "FIG_F6_GAIN_CHANGES_BY_BUS.png"):
        assert (HERE / name).stat().st_size > 10000
    report = (HERE / "FINAL_REPORT.md").read_text()
    assert f"FINAL_BEST_VALIDATED_TAU_CRIT_MS: {float(best['tau_crit_ms']):.9f}" in report
    assert "OPTIMIZATION_CONVERGED: NO" in report
    assert "POSTER_READY: NO" in report
    print(f"DELIVERY_CHECK_PASS best={cid} tau_ms={float(best['tau_crit_ms']):.9f} validated_designs={len(validated)}")


if __name__ == "__main__":
    main()
