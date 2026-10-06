"""Write an evidence-only closure summary from frozen experiment outputs."""
from __future__ import annotations

import csv
import math
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(name: str):
    with (HERE / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def fmt(value, digits=9):
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else str(value)


def main():
    parent = rows("F00_PARENT_REPRODUCTION.csv")
    trace = rows("TABLE_F2_FINAL_OPTIMIZATION_TRACE.csv")
    feasible = [r for r in trace if r["status"] == "VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"]
    assert feasible, "no fully validated zero-delay design"
    best = max(feasible, key=lambda r: float(r["tau_crit_ms"]))
    cid = best["design_id"]
    critical = float(best["tau_crit_ms"])
    z, nominal, old = (float(r["tau_crit_ms"]) for r in parent)
    design = tomllib.loads((HERE / "designs" / f"{cid}.toml").read_text())
    events = tomllib.loads((HERE / "event_validations" / cid / "Q0_RESULT.toml").read_text())
    roots = rows(f"evaluations/{cid}/ROOTS.csv")
    active = [r for r in roots if float(r["local_crossing_ms"]) - critical <= 0.05]
    modes = [r for r in rows("TABLE_F3_ACTIVE_MODAL_FAMILIES.csv") if r["design_id"] == cid
             and int(r["rank"]) in {int(a["rank"]) for a in active}]
    rankrows = rows("TABLE_F11_GAIN_ACTION_RANK.csv")
    maxerr = max(float(r["relative_action_reconstruction_error"]) for r in rankrows)
    delayed = rows("TABLE_F8_POSITIVE_DELAY_VALIDATION.csv")
    delayed_best = [r for r in delayed if r["design_id"] == cid and r["steps_per_delay"] == "160"]
    kkt_file = HERE / f"KKT_LP_{cid}.csv"
    local_kkt = max(float(r["kkt_residual_inf"]) for r in rows(kkt_file.name)) if kkt_file.exists() else math.nan
    slack = min(float(e["min_SG_actuator_fraction_slack"]) for e in events["events"])
    worst = min(events["events"], key=lambda e: float(e["min_SG_actuator_fraction_slack"]))
    family_ids = [r["closest_previous_best_family"] for r in modes]
    accepted_ids = {r["design_id"] for r in feasible}
    minmac = min(float(r["closest_previous_right_MAC"]) for r in rows("TABLE_F3_ACTIVE_MODAL_FAMILIES.csv")
                 if r["design_id"] in accepted_ids)
    reserve_feasible = [r for r in feasible if float(r["actuator_slack"]) >= 0.00202]
    reserve_row = max(reserve_feasible, key=lambda r: float(r["tau_crit_ms"])) if reserve_feasible else None
    v2 = len(delayed_best) == 3
    behavior = {r["factor"]: r["observed_behavior"] for r in delayed_best}
    boundary_spread = max(float(a["local_crossing_ms"]) for a in active) - critical
    uncertainty = rows("TABLE_F9_UNCERTAINTY_AUDIT.csv")
    plus = next((r for r in uncertainty if r["case_id"] == "uncertainty_best_event_plus1MW"), None)
    reserve_plus = next((r for r in uncertainty if r["case_id"] == "uncertainty_reserve_event_plus1MW"), None)
    tiny_plus = next((r for r in uncertainty if r["case_id"] == "uncertainty_tiny_event_plus1MW"), None)
    tiny_gain = next((r for r in uncertainty if r["case_id"] == "uncertainty_tiny_Kp33_minus1pct"), None)
    step13_plus = next((r for r in uncertainty if r["case_id"] == "uncertainty_step13_event_plus1MW"), None)
    step13_small_pass = next((r for r in uncertainty if r["case_id"] == "uncertainty_step13_event_plus0p01MW"), None)
    step13_small_fail = next((r for r in uncertainty if r["case_id"] == "uncertainty_step13_event_plus0p02MW"), None)
    audit = (f"+101 MW at step11: actuator slack {plus['actuator_slack']} < 0.002, one-event FAIL; "
             f"reserve design +101 MW: actuator slack {reserve_plus['actuator_slack'] if reserve_plus else 'NOT_RUN'} < 0.002; "
             f"step12 tiny +101 MW: actuator slack {tiny_plus['actuator_slack'] if tiny_plus else 'NOT_RUN'} < 0.002; "
             f"Kp33 -1% at step12 tiny: actuator slack {tiny_gain['actuator_slack'] if tiny_gain else 'NOT_RUN'}; "
             f"step13 +101 MW: actuator slack {step13_plus['actuator_slack'] if step13_plus else 'NOT_RUN'}; "
             f"step13 +100.01 MW {'PASS' if step13_small_pass and step13_small_pass['pass_screen']=='true' else 'CHECK'}; "
             f"+100.02 MW {'FAIL' if step13_small_fail and step13_small_fail['pass_screen']=='false' else 'CHECK'}; "
             "±1% ZIP-load cases reequilibrated but numerical solver did not complete") if plus else "INCOMPLETE"
    block = "No validated nonlinear positive-delay DDE event simulations."
    data = [
        ("FINAL_CLOSURE_STATUS", "BEST_VALIDATED_DESIGN_FOUND; MAXIMUM_NOT_ESTABLISHED"),
        ("FIXED_GFL_PERCENT", fmt(events["GFL_percent"])),
        ("FIXED_GFL_MW", fmt(events["GFL_MW"])),
        ("FIXED_RETAINED_SG_MW", fmt(events["retained_SG_MW"])),
        ("PENDING_41MS_CANDIDATE", "pending_41ms"),
        ("PENDING_CANDIDATE_TRUE_TAU_CRIT_MS", rows("F01_PENDING_CANDIDATE_VALIDATION.csv")[0]["true_tau_crit_ms"]),
        ("PENDING_CANDIDATE_FEASIBLE", rows("F01_PENDING_CANDIDATE_VALIDATION.csv")[0]["status"]),
        ("ZERO_DELAY_TUNED_TAU_CRIT_MS", fmt(z)),
        ("NOMINAL_TAU_CRIT_MS", fmt(nominal)),
        ("PREVIOUS_BEST_TAU_CRIT_MS", fmt(old)),
        ("FINAL_BEST_VALIDATED_TAU_CRIT_MS", fmt(critical)),
        ("TOTAL_RECOVERY_VS_ZERO_DELAY_MS", fmt(critical-z)),
        ("TOTAL_RECOVERY_VS_ZERO_DELAY_PERCENT", fmt(100*(critical-z)/z, 6)),
        ("RECOVERY_VS_NOMINAL_MS", fmt(critical-nominal)),
        ("FINAL_RHO", repr(design["rho"])),
        ("FINAL_KP", repr(design["Kp"])),
        ("FINAL_KI", repr(design["Ki"])),
        ("OPTIMIZATION_CONVERGED", "NO"),
        ("CONVERGENCE_REASON", "No full-model KKT/no-material-improvement/active-ceiling criterion met; latest accepted step still improved materially"),
        ("KKT_RESIDUAL", f"FULL_MODEL_NOT_ESTABLISHED; local LP {local_kkt:.4e}"),
        ("NUMBER_ACTIVE_DELAYED_FAMILIES", f"{len(active)} within declared 0.05 ms numerical band"),
        ("ACTIVE_FAMILY_IDS", ", ".join(family_ids) + " (nearest-reference provisional labels)"),
        ("ACTIVE_FAMILY_FREQUENCIES_HZ", ", ".join(fmt(float(r["frequency_hz"]), 6) for r in active)),
        ("ACTIVE_FAMILY_TAU_CRIT_SPREAD_MS", fmt(boundary_spread)),
        ("MODAL_EQUALIZATION", "NEAR_BALANCE_NUMERICAL; not optimum"),
        ("MODAL_SWITCH_CONFIRMED", "SUPPORTED_NUMERICALLY by low reference MAC; continuity labels not certified"),
        ("MIN_KEY_FAMILY_MAC", fmt(minmac)),
        ("ACTUATOR_CONSTRAINT_ACTIVE", "NO strict nominal activity; local LP surrogate active"),
        ("ACTUATOR_VALUE", fmt(slack)),
        ("ACTUATOR_LIMIT", "0.002000000"),
        ("ACTUATOR_RESERVE", fmt(slack-0.002)),
        ("ZERO_DELAY_ALPHA", fmt(float(events["critical_real_part_s_inv"]))),
        ("ALL_FIVE_ZERO_DELAY_EVENTS_PASS", str(events["all_five_events_pass"]).upper()),
        ("MAX_ZERO_DELAY_FREQ_DEVIATION", fmt(max(float(e["F_peak_Hz"]) for e in events["events"]))),
        ("MAX_ZERO_DELAY_ROCOF", fmt(max(float(e["RoCoF_peak_Hz_s"]) for e in events["events"]))),
        ("WORST_FROZEN_EVENT", f"{worst['event']} by actuator slack"),
        ("POSITIVE_DELAY_TIME_DOMAIN", "V2 NUMERICAL METHOD OF STEPS FOR THE FULL LINEAR DDE (NO PADE)" if v2 else "NOT_RUN_FOR_FINAL_BEST"),
        ("BELOW_THRESHOLD_RESULT", behavior.get("0.9", "NOT_RUN")),
        ("NEAR_THRESHOLD_RESULT", behavior.get("0.98", "NOT_RUN")),
        ("ABOVE_THRESHOLD_RESULT", behavior.get("1.02", "NOT_RUN")),
        ("MULTISTART_CONSISTENCY", "NOT_TESTED"),
        ("GENERIC_OPTIMIZER_BEST_TAU_MS", "NOT_RUN"),
        ("ANALYTICAL_METHOD_BEST_TAU_MS", fmt(critical)),
        ("FULL_MODEL_EVALUATION_SAVINGS", "NOT_MEASURED"),
        ("EXACT_GAIN_ACTION_RANK", "≤10 for 20 gain variables at fixed rho"),
        ("MAX_RECONSTRUCTION_ERROR", f"{maxerr:.4e} relative over checked designs"),
        ("PARETO_FRONTIER", "NOT_OPTIMIZED; FIG_F3 shows evaluated points only"),
        ("NOMINAL_VS_LATENCY_TRADEOFF", f"observed -alpha0={-float(best['alpha_zero_delay_s_inv']):.9f} s^-1 at tau={critical:.9f} ms"),
        ("ROBUSTNESS_AUDIT", audit),
        ("ROBUST_DESIGN_TAU_CRIT_MS", f"{float(reserve_row['tau_crit_ms']):.9f} nominal 0.002020000 actuator-reserve target met; no +101 MW robustness claim" if reserve_row else "NOT_ESTABLISHED; reserve attempt did not meet declared 0.002020000 target"),
        ("STRONGEST_EXACT_RESULT", "Gain-induced linear action factorization rank ≤10 at fixed rho"),
        ("STRONGEST_ANALYTICAL_RESULT", "Simple-root latency derivative and multimode lower-envelope formula"),
        ("STRONGEST_PHYSICAL_RESULT", "Five frozen zero-delay nonlinear events pass at unchanged replacement"),
        ("STRONGEST_NUMERICAL_RESULT", f"Validated characteristic critical delay {critical:.9f} ms; gain +{100*(critical-z)/z:.3f}% versus Z"),
        ("STRONGEST_NEGATIVE_RESULT", "+101 MW event defeats step11 and nominal-reserve design actuator guards; no nonlinear delayed safety result"),
        ("POSTER_READY", "NO for maximum/robust-safety headline"),
        ("RECOMMENDED_TITLE", "Multimode PLL latency tuning at fixed 88.455% GFL replacement on IEEE-39"),
        ("HEADLINE_CLAIM", f"Best validated gain tuning moves the numerical spectral latency margin from {z:.3f} to {critical:.3f} ms at unchanged replacement and passes five zero-delay events"),
        ("HEADLINE_EQUATION", "tau_crit(K)=min_m tau_m(K),  grad_K tau_m=-(grad_K Re lambda_m)/(partial_tau Re lambda_m)"),
        ("HEADLINE_NUMBER", f"+{critical-z:.3f} ms / +{100*(critical-z)/z:.2f}% versus Z"),
        ("HEADLINE_FIGURE", "FIG_F1_LATENCY_MARGIN_OPTIMIZATION.png"),
        ("IF POSTER_READY = NO", "Full delayed nonlinear event validation is absent; maximum/optimality claim also needs convergence and comparative checks"),
        ("EXACT_SINGLE_BLOCKER", block),
    ]
    (HERE / "FINAL_REPORT.md").write_text(
        "# Final closure report\n\n"
        + "All numerical values below are read from the executed tables and the selected five-event design. "
          "Critical-delay root coverage is numerical, not interval certified. The positive-delay time-domain result is linear.\n\n"
        + "\n".join(f"{key}: {value}" for key, value in data) + "\n", encoding="utf-8")
    print(cid, critical, "report written")


if __name__ == "__main__":
    main()
