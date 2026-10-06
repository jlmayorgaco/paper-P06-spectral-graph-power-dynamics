"""Write the claim ledger only from executed, checked experiment artifacts."""

from __future__ import annotations

import csv
import json
import math
import shutil
import statistics
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent


def rows(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def result(cid: str) -> dict:
    return tomllib.loads((HERE / "event_validations" / cid / "Q0_RESULT.toml").read_text())


def design(cid: str) -> dict:
    path = HERE / "designs" / ("Z_zero_delay_tuned.toml" if cid == "Z" else
                             "N_nominal.toml" if cid == "N" else f"{cid}.toml")
    return tomllib.loads(path.read_text())


def main() -> None:
    trace = rows("L3_OPTIMIZATION_TRACE.csv")
    accepted = [r for r in trace if r["accepted"] == "True"]
    assert accepted, "no fully validated design"
    best = max(accepted, key=lambda r: float(r["tau_crit_ms"]))
    cid = best["candidate_id"]
    best_event = result(cid)
    assert best_event["all_five_events_pass"] and len(best_event["events"]) == 5
    assert best["spectral_bracket_pass"] == "True"
    d = design(cid)
    n = next(r for r in trace if r["candidate_id"] == "N")
    z = next(r for r in trace if r["candidate_id"] == "Z")
    increment_n = float(best["tau_crit_ms"]) - float(n["tau_crit_ms"])
    increment_z = float(best["tau_crit_ms"]) - float(z["tau_crit_ms"])
    l1 = rows("L1_ACTION_SPACE_VALIDATION.csv")
    l2 = rows("L2_TAU_MARGIN_GRADIENT_VALIDATION.csv")
    sw = rows("L4_MODE_SWITCH_POINT.csv")[0]
    repl = rows("L0_REPLACEMENT_ONLY_REPRODUCTION.csv")
    repl_effect = float(repl[1]["tau_crit_ms"]) - float(repl[0]["tau_crit_ms"])
    worst_frequency = max(best_event["events"], key=lambda e: float(e["F_peak_Hz"]))
    worst_rocof = max(best_event["events"], key=lambda e: float(e["RoCoF_peak_Hz_s"]))
    limiting_act = min(best_event["events"], key=lambda e: float(e["min_SG_actuator_fraction_slack"]))
    coverage = HERE / "evaluations" / cid / "ROOT_COVERAGE.csv"
    roots = rows(f"evaluations/{cid}/ROOTS.csv")
    crossing_values = sorted(float(r["local_crossing_ms"]) for r in roots)
    counts = rows(f"evaluations/{cid}/ROOT_COUNTS.csv")
    unsafe_count = max(int(r["root_count"]) for r in counts if r["status"] == "UNSAFE")
    covered_count = 2 * sum(float(r["local_crossing_ms"]) <= max(float(c["tau_ms"]) for c in counts if c["status"] == "UNSAFE") for r in roots)
    root_cover_status = (rows(f"evaluations/{cid}/ROOT_COVERAGE.csv")[0]["status"] if coverage.exists()
                         else "INITIAL_DISCOVERY_MATCHES_CONTOUR" if unsafe_count == covered_count else "INCOMPLETE")
    assert root_cover_status != "INCOMPLETE", (cid, unsafe_count, covered_count)
    next_prediction_path = HERE / f"L3_FOLLOWUP_PREDICTION_{cid}.csv"
    next_prediction = (rows(next_prediction_path.name)[0] if next_prediction_path.exists() else None)
    shutil.copy2(HERE / "designs" / f"{cid}.toml", HERE / "BEST_FOUND_DESIGN.toml")
    summary = {
        "claim_status": "FULL_MODEL_NUMERICAL_RESULT_NOT_GLOBAL_OPTIMUM",
        "fixed_GFL_percent": float(best_event["GFL_percent"]),
        "fixed_GFL_MW": float(best_event["GFL_MW"]),
        "fixed_retained_SG_MW": float(best_event["retained_SG_MW"]),
        "best_fully_validated_found_id": cid,
        "best_local_numerical_uniform_latency_threshold_ms": float(best["tau_crit_ms"]),
        "nominal_N_threshold_ms": float(n["tau_crit_ms"]),
        "zero_delay_tuned_Z_threshold_ms": float(z["tau_crit_ms"]),
        "improvement_vs_N_ms": increment_n,
        "improvement_vs_Z_ms": increment_z,
        "replacement_only_fixed_nominal_gain_effect_ms": repl_effect,
        "best_alpha_zero_delay_s_inv": float(best["alpha_zero_delay"]),
        "best_critical_frequency_hz": float(best["critical_delay_frequency_hz"]),
        "best_min_SG_actuator_slack": float(limiting_act["min_SG_actuator_fraction_slack"]),
        "limiting_actuator_event": limiting_act["event"],
        "best_max_frequency_deviation_Hz": float(worst_frequency["F_peak_Hz"]),
        "best_max_RoCoF_Hz_s": float(worst_rocof["RoCoF_peak_Hz_s"]),
        "best_Kp": d["Kp"], "best_Ki": d["Ki"], "rho": d["rho"],
        "gain_action_rank": 10,
        "gain_action_validation_evaluations": len(l1),
        "max_relative_action_reconstruction_error": max(float(r["relative_reconstruction_error"]) for r in l1),
        "baseline_gradient_sign_agreements": sum(r["sign_agreement"] == "true" for r in l2),
        "baseline_gradient_tests": len(l2),
        "baseline_gradient_median_relative_error": statistics.median(float(r["relative_error"]) for r in l2),
        "baseline_gradient_max_relative_error": max(float(r["relative_error"]) for r in l2),
        "modal_switch_eta": float(sw["eta_estimate"]),
        "modal_switch_right_MAC": float(sw["physical_right_vector_MAC"]),
        "root_coverage_status": root_cover_status,
        "best_first_three_mode_crossing_gaps_ms": [x-crossing_values[0] for x in crossing_values[:3]],
        "delayed_nonlinear_validation": "NOT_RUN",
        "heterogeneous_delay_optimization": "NOT_RUN",
        "global_optimality_bound": "NOT_AVAILABLE",
        "local_optimizer_stationarity": "NOT_ESTABLISHED_SEARCH_STOPPED_BEFORE_CONVERGENCE",
        "next_unvalidated_local_predictor_ms": (float(next_prediction["predicted_tau_ms"])
                                                 if next_prediction else None),
        "numerical_DDE_root_count_interval_certificate": "NOT_AVAILABLE",
    }
    (HERE / "RESULT_SUMMARY.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")

    (HERE / "README.md").write_text(f"""# Latency-robust PLL co-design at fixed IEEE-39 replacement

This self-contained **experiment output** reuses the frozen repository model; see `PROVENANCE.md` and `EXPERIMENT_MANIFEST.json`. At fixed {summary['fixed_GFL_percent']:.6f}% GFL replacement, it tunes ten independent PLL proportional and integral gains against a common, exogenous delay on each PLL measurement/error path. No files from prior experiments were edited.

The best fully validated design **found in the executed search** is `{cid}`. Its numerical local uniform-delay threshold is **{summary['best_local_numerical_uniform_latency_threshold_ms']:.6f} ms**, compared with **{summary['nominal_N_threshold_ms']:.6f} ms** for nominal N and **{summary['zero_delay_tuned_Z_threshold_ms']:.6f} ms** for zero-delay-tuned Z. It passes the full zero-delay physical spectrum and all five frozen nonlinear 61-second events. The corresponding `rho`, `Kp`, and `Ki` vectors are in `BEST_FOUND_DESIGN.toml`.

Read `STATUS.md` for claim gates, `THEORY_LATENCY_ROBUST_CODESIGN.md` for the exact rank-10 action and simple-root sensitivity, `L3_OPTIMIZATION_TRACE.csv` for accepted and rejected steps, `L4_MODAL_FAMILY_MAP.csv` for modal identity, and `POSTER_CLAIM_LEDGER.md` for the allowed claims. `figures/` contains diagnostic plots generated from the CSVs.

The word “exact” describes the exponential DDE characteristic and algebraic gain action. Numerical roots, contour counts, and optimization are Float64 calculations. The threshold is a **numerically observed local first crossing**, without an interval certificate of all characteristic roots or a global optimality proof. The nonlinear events use **zero delay**; the experiment does not establish nonlinear safety under positive delay.
""", encoding="utf-8")

    (HERE / "STATUS.md").write_text(f"""# Experiment status and gates

| Gate | Result | Evidence |
|---|---|---|
| L0 baseline parity | PASS | `L0_REPRODUCTION.csv`; `L0_REPLACEMENT_ONLY_REPRODUCTION.csv` |
| L1 fixed-ρ gain action | PASS (exact identity, Float64 check) | `L1_ACTION_SPACE_VALIDATION.csv`: rank 10 in {len(l1)} evaluations; max relative error {summary['max_relative_action_reconstruction_error']:.3g} |
| L2 simple-root gradient | PASS for tested smooth cases | {summary['baseline_gradient_sign_agreements']}/{summary['baseline_gradient_tests']} signs, median relative error {summary['baseline_gradient_median_relative_error']:.3g}; near-switch envelope nonsmooth in `L2_INTERMEDIATE_GRADIENT_VALIDATION.csv` |
| L3 two-start co-design | EXECUTED, best found; not converged | `L3_OPTIMIZATION_TRACE.csv`; `{cid}` passes five events; prescribed <0.01 ms/KKT stopping gate not reached |
| L4 family switch | SUPPORTED NUMERICALLY | `L4_MODE_SWITCH_POINT.csv`; η≈{summary['modal_switch_eta']:.6f}, MAC≈{summary['modal_switch_right_MAC']:.4f} |
| Complete numerical contour around best threshold | PASS | `evaluations/{cid}/ROOT_COUNTS.csv`; coverage: {root_cover_status} |
| Positive-delay nonlinear DDE events | NOT RUN | No delayed nonlinear simulator in this experiment |
| Global optimum/interval root certificate | NOT AVAILABLE | No validated global bound or interval arithmetic |

The frequency and actuator event limits were frozen in `EXPERIMENT_MANIFEST.json`. The rejected full second step raised the spectral threshold but failed the bus16 +100 MW actuator guard (`event_screens/N_step1_multimode_lp_next_full/`). The winning design has only {summary['best_min_SG_actuator_slack']-0.002:.3g} actuator-fraction slack above that guard, so it is sensitive to numerical/model uncertainty.
The next frozen-trust LP predictor, if present, is **unvalidated** and is excluded from the best-found result. Its predicted threshold is {summary['next_unvalidated_local_predictor_ms']} ms; see `PENDING_CANDIDATES.md`.
""", encoding="utf-8")

    (HERE / "PENDING_CANDIDATES.md").write_text(f"""# Unvalidated continuation

After the last fully validated design `{cid}`, the same local analytical LP produced a new candidate with predicted threshold {summary['next_unvalidated_local_predictor_ms']} ms. It has **no exact spectral corrector and no frozen-event validation** in this archive. It is not included in the feasible result or poster claims. This material predicted gain is direct evidence that the local optimizer has not met its stopping criterion. The candidate TOML and predictor CSV are preserved only to resume the search.
""", encoding="utf-8")

    (HERE / "POSTER_CLAIM_LEDGER.md").write_text(f"""# Poster claim ledger — executed evidence only

1. **EXACT_IDENTITY:** At fixed ρ, the 20 PLL gain variables act through ten measurement channels, so the descriptor gain update has rank ≤10. Forty-eight full-matrix checks found rank 10 and maximum relative reconstruction error {summary['max_relative_action_reconstruction_error']:.3g}. See `THEORY_LATENCY_ROBUST_CODESIGN.md` and `L1_ACTION_SPACE_VALIDATION.csv`.
2. **NUMERICALLY_VALIDATED:** The simple-root analytical latency derivative agreed in sign with all {len(l2)} tested finite differences (median relative error {summary['baseline_gradient_median_relative_error']:.3g}). It is not a derivative of the nonsmooth two-family envelope at a switch. See L2 tables.
3. **SUPPORTED_EMPIRICAL:** At the same {summary['fixed_GFL_percent']:.6f}% GFL replacement, Z and N had thresholds {summary['zero_delay_tuned_Z_threshold_ms']:.3f} and {summary['nominal_N_threshold_ms']:.3f} ms. The two limiting physical right eigenvectors have MAC {float(rows('L0_MODAL_OVERLAP.csv')[0]['full_physical_right_vector_MAC']):.4f}. The tracked families co-limit at η≈{summary['modal_switch_eta']:.6f} along the frozen Z→N log-gain path. See L0/L4 tables and `figures/FIG_L4_MODAL_FAMILY_SWITCH.png`.
4. **FULL_MODEL_NUMERICAL_RESULT:** The best **found**, fully zero-delay-event-validated gain set has numerical local uniform-delay threshold {summary['best_local_numerical_uniform_latency_threshold_ms']:.3f} ms, an improvement of {increment_n:.3f} ms over N and {increment_z:.3f} ms over Z, at unchanged ρ and {summary['fixed_retained_SG_MW']:.2f} MW retained SG. All five frozen zero-delay events pass. See `BEST_FOUND_DESIGN.toml`, `RESULT_SUMMARY.json`, and `L3_OPTIMIZATION_TRACE.csv`.
5. **NEGATIVE_RESULT:** On the independently re-evaluated fixed-nominal-gain 87.5→{summary['fixed_GFL_percent']:.3f}% replacement path, the observed latency change is only {repl_effect:+.4f} ms. This corridor does not support a material replacement-only latency effect.
6. **NEGATIVE_RESULT:** A more aggressive predictor passed the spectral check but failed the bus16 +100 MW SG-actuator constraint. A spectral-only gain is not physically accepted. See `event_screens/N_step1_multimode_lp_next_full/`.
7. **NUMERICALLY_VALIDATED METHOD CHECK:** The first root seed can miss the earliest of three fast pairs. The full contour and a finer frequency scan corrected later candidate thresholds before acceptance; both initial and corrected tables are retained in `evaluations/`. Never infer the margin from a single tracked root.
8. **SUPPORTED_LOCAL:** At the best found design, the next two discovered delayed-mode crossings lie {crossing_values[1]-crossing_values[0]:.6f} and {crossing_values[2]-crossing_values[0]:.6f} ms above the first. This is local numerical co-criticality, not a proof of a unique optimum. See `evaluations/{cid}/ROOTS.csv` and `figures/FIG_L4_COCRITICAL_MODE_GAPS.png`.

**Do not claim:** global optimality; interval-certified exact DDE spectrum; heterogeneous-delay optimality; positive-delay nonlinear disturbance safety; or robustness to uncertainty in the near-active SG actuator guard. These are untested.
""", encoding="utf-8")

    (HERE / "REPRODUCE.md").write_text("""# Reproduction

For an independent **best-design recheck**, run from the repository root with the Julia project and Python packages recorded in `PROVENANCE.md`:

```powershell
julia --project=. experiments/latency_robust_pll_codesign_20261003/evaluate_design.jl experiments/latency_robust_pll_codesign_20261003/BEST_FOUND_DESIGN.toml
julia --project=. experiments/latency_robust_pll_codesign_20261003/run_event_local.jl experiments/latency_robust_pll_codesign_20261003/BEST_FOUND_DESIGN.toml
```

For the entire campaign, use the archived stage scripts in order: L0 reproduction and both event validations; L1 frozen probes and action validation; L2 gradients, finite differences and family scan; L3 initial predictor/corrector and each design's `evaluate_design.jl`, root-coverage discovery, and five-event correction; L4 family map; then `finalize_results.py` and `write_report.py`. `L3_OPTIMIZATION_TRACE.csv` and each `designs/*.toml` give the executed candidate order. Some freeze scripts refuse to overwrite existing artifacts: rerun the full pipeline in a **separate clean copy**, retaining `EXPERIMENT_MANIFEST.json` and the frozen probe tables. Do not regenerate preregistration after seeing results.

The original repository structure and prior experiment modules are required. `source_snapshot/` and `input_snapshot/` preserve source and data inputs captured after the runs for comparison or recovery in a separate checkout. A five-event validation takes several minutes. Rechecking with `BEST_FOUND_DESIGN.toml` writes a separate `BEST_FOUND_DESIGN` output directory under this experiment.
""", encoding="utf-8")

    (HERE / "FINAL_REPORT.md").write_text(f"""# Final experiment report

**Question.** At a fixed {summary['fixed_GFL_percent']:.6f}% GFL replacement in IEEE-39, can independent PLL gains enlarge the numerical uniform measurement-latency margin while preserving the frozen zero-delay dynamics and five external events?

**Observed answer.** Yes, locally. The best fully validated gain set found has τ≈{summary['best_local_numerical_uniform_latency_threshold_ms']:.6f} ms, versus N={summary['nominal_N_threshold_ms']:.6f} ms and Z={summary['zero_delay_tuned_Z_threshold_ms']:.6f} ms. The gain vectors are in `BEST_FOUND_DESIGN.toml`; ρ stays fixed, with {summary['fixed_GFL_MW']:.2f} MW GFL and {summary['fixed_retained_SG_MW']:.2f} MW SG. This is an executed feasible improvement, not a proof of the maximum over all gains.

**Mechanism.** Twenty gains update the characteristic matrix through ten PLL measurement channels exactly. The simple-root sensitivity accurately predicts tested local directions, but the margin envelope switches between two nearly orthogonal physical modal families at η≈{summary['modal_switch_eta']:.6f}; the root MAC there is {summary['modal_switch_right_MAC']:.4f}. The exact exponential characteristic and numerical contour are the spectral corrector. The strongest limit on the found tuning is the zero-delay SG actuator under `{limiting_act['event']}`: fraction slack {float(limiting_act['min_SG_actuator_fraction_slack']):.9f} against a 0.002 floor.

**Frozen event evidence.** All five 61-second zero-delay events pass. Maximum |Δf| is {float(worst_frequency['F_peak_Hz']):.6f} Hz (event `{worst_frequency['event']}`); maximum RoCoF is {float(worst_rocof['RoCoF_peak_Hz_s']):.6f} Hz/s (event `{worst_rocof['event']}`). Full event-level voltage, DC, current-ratio, and actuator observations are stored in `event_validations/{cid}/Q0_RESULT.toml`. A current ratio is reported, but no independent current-limit safety claim is made.

**Falsification.** Holding gains nominal, the independently recomputed 87.5%→{summary['fixed_GFL_percent']:.3f}% replacement path changes the observed threshold by only {repl_effect:+.6f} ms. A larger gain step reached a higher spectral threshold but breached the frozen SG-actuator floor. Near the family switch the lower-envelope gradient is nonsmooth, so a single active-mode derivative is not a valid global direction. These results constrain the poster narrative.

**Validity.** The matrix update and determinant lemma are identities in the declared model. Spectrum/contour, event simulation, gradients, and modal labels are Float64 numerical evidence. The reported threshold is local and observed over the checked delay range, without exhaustive DDE root certification. The local search still had material improving steps and did not meet the prescribed <0.01 ms/KKT stopping gate. We did not execute nonlinear positive-delay DDE events, heterogeneous-delay optimization, uncertainty sweeps, a KKT stationarity proof, or a global upper bound. Therefore neither delayed dynamic security nor globally optimal tuning is established.

**Poster decision.** The rank-10 action, measured gradient accuracy, and two-family switch form a concrete mathematical result. A strong latency-robust operating claim requires positive-delay nonlinear events and a wider robustness check at the nearly active actuator limit. The poster should say “best found” rather than “optimal.”
""", encoding="utf-8")
    print("REPORT_DONE", cid, summary["best_local_numerical_uniform_latency_threshold_ms"])


if __name__ == "__main__":
    main()
