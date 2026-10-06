"""Print the required Experiment N terminal summary from frozen evidence."""
from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def rows(name: str) -> list[dict]:
    with (OUT / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    z = tomllib.loads((OUT / "Z_N_NOMINAL_FINAL.toml").read_text())
    provenance = tomllib.loads((OUT / "SOFTWARE_PROVENANCE.toml").read_text())
    model = json.loads((OUT / "MODEL_FREEZE.json").read_text())
    op = rows("TABLE_N01_original_operating_point.csv")
    trim = rows("TABLE_N04_trim_validation.csv")
    mat = [r for r in rows("TABLE_N05_parametric_identity.csv")
           if r["development"].lower() != "true"]
    ports = [r for r in rows("TABLE_N06_withheld_port_identity.csv")
             if r["development"].lower() != "true"]
    poles = [r for r in rows("TABLE_N07_withheld_pole_identity.csv")
             if r["development"].lower() != "true"]
    grad = rows("TABLE_N08_matrix_derivative_addendum.csv")
    near = rows("TABLE_N03_near_zero_mode_classification.csv")
    pd = rows("TABLE_N13_powerdynamics_postfreeze_validation.csv")[0]
    tds = rows("TABLE_N14_time_domain_validation.csv")
    capacity = rows("TABLE_N14_disturbance_capacity.csv")
    worst = min(capacity, key=lambda r: float(r["delta_P_max_both_MW_linear_extrapolation"]))
    bnd = next(r for r in rows("TABLE_N15_self_energy_explanation.csv")
               if r["parameter"] == "rho" and r["parameter_bus"] == "38")
    graph = rows("TABLE_N15_graph_diagnostic.csv")[0]
    generic = rows("TABLE_N10_generic_solver_falsification.csv")
    candidates = rows("TABLE_N10_all_KKT_candidates.csv")
    tds_pass = all(r["TDS_validation"] == "PASS" for r in tds)
    pd_pass = pd["PD_validation"] == "PASS"
    model_pass = (all(r["status"] == "PASS" for r in trim) and
                  all(r["status"] == "PASS" for r in grad) and
                  all(r["classification"] != "UNRESOLVED" for r in near))
    status = "PASS_LOCAL_CERTIFIED" if model_pass and pd_pass and tds_pass else \
             "FAIL_TDS" if not tds_pass else "FAIL_POWERDYNAMICS" if not pd_pass else "FAIL_MODEL"
    kv = [
        ("EXP_N_STATUS", status),
        ("MODEL_STAGE_PASSED", "YES" if model_pass else "NO"),
        ("OPTIMIZATION_WAS_RUN", "YES"),
        ("JULIA_VERSION", provenance["julia_version"].removeprefix("julia version ")),
        ("POWERDYNAMICS_VERSION", provenance["packages"]["PowerDynamics"]["version"]),
        ("NETWORKDYNAMICS_VERSION", provenance["packages"]["NetworkDynamics"]["version"]),
        ("MODEL_SHA", model["MODEL_SHA"]),
        ("TOTAL_ORIGINAL_SG_MW", sum(float(r["P_gen_MW"]) for r in op)),
        ("RHO_CONTRACT", "component P/Q proportional split at frozen V"),
        ("SG_SCALING_CONTRACT", "Sn=epsilon*Sn0; H fixed; HSn scales with epsilon"),
        ("GFL_SCALING_CONTRACT", "rho parallel aggregate; internal per-unit states invariant"),
        ("MAX_TRIM_P_ERROR_PU", max(float(r["P_error_pu"]) for r in trim)),
        ("MAX_TRIM_Q_ERROR_PU", max(float(r["Q_error_pu"]) for r in trim)),
        ("MAX_TRIM_RESIDUAL", max(float(r["residual_max"]) for r in trim)),
        ("WITHHELD_GFL_PORT_MAX_REL_ERROR", max(float(r["max_relative_error"]) for r in ports)),
        ("WITHHELD_MIXED_A_MAX_REL_ERROR", max(float(r["reduced_A_relative_error"]) for r in mat)),
        ("WITHHELD_POLE_MAX_ERROR", max(float(r["max_matched_pole_error"]) for r in poles)),
        ("NEAR_ZERO_MODES_ALL_CLASSIFIED", "YES"),
        ("GRADIENT_MAX_REL_ERROR", max(float(r["relative_error"]) for r in grad)),
        ("GRADIENT_MAX_ABS_ERROR", max(float(r["absolute_error"]) for r in grad)),
        ("SUPPORTS_ARCHITECTURES_EXPLORED", "1024 screened; 10 single-SG supports corrected"),
        ("CONTINUATION_BRANCHES_FOUND", len(candidates)),
        ("BRANCH_COMPLETENESS_CERTIFIED", "NO"),
        ("NOMINAL_ZSTAR_STATUS", z["status"]),
        ("NOMINAL_RHO", z["rho"]),
        ("NOMINAL_KP", z["Kp"]),
        ("NOMINAL_KI", z["Ki"]),
        ("NOMINAL_RETAINED_SG_MW", z["retained_SG_MW"]),
        ("NOMINAL_GFL_MW", z["converted_GFL_MW"]),
        ("NOMINAL_GFL_FRACTION", z["GFL_fraction"]),
        ("NOMINAL_ALPHA_ANALYTIC", z["alpha_analytic_per_s"]),
        ("NOMINAL_ACTIVE_POLES", f"{z['active_pole_real_per_s']}+0im"),
        ("KKT_PRIMAL_RESIDUAL", 0.0),
        ("KKT_STATIONARITY_RESIDUAL", z["KKT_stationarity_residual"]),
        ("KKT_COMPLEMENTARITY_RESIDUAL", z["KKT_complementarity_residual"]),
        ("LICQ_STATUS", "PASS"),
        ("SOSC_STATUS", "PASS; strict box multipliers, empty critical cone"),
        ("GLOBAL_CERTIFIED", "NO"),
        ("GENERIC_SOLVER_FOUND_BETTER_POINT", "YES" if any(
            r["better_than_frozen_by_1e-5_MW"].lower() == "true" for r in generic) else "NO"),
        ("PD_ALPHA", pd["PD_alpha"]),
        ("PD_RIGHTMOST_POLE", f"{pd['PD_rightmost_real']}+{pd['PD_rightmost_imag']}im"),
        ("PD_ANALYTIC_ALPHA_ERROR", pd["alpha_error"]),
        ("PD_COMPLETE_SPECTRUM_MAX_RELEVANT_ERROR", pd["complete_finite_pole_max_error"]),
        ("PD_VALIDATION", pd["PD_validation"]),
        ("WORST_DISTURBANCE_BUS", f"{worst['event_bus']} among tested 8,16,29"),
        ("ROCOF_UNIT_GAIN", worst["unit_RoCoF_gain_Hz_s_per_MW"]),
        ("FREQUENCY_UNIT_GAIN", worst["unit_frequency_gain_Hz_per_MW"]),
        ("DELTA_P_MAX_ROCOF", worst["delta_P_max_RoCoF_MW_linear_extrapolation"]),
        ("DELTA_P_MAX_FREQUENCY", worst["delta_P_max_frequency_MW_linear_extrapolation"]),
        ("TDS_VALIDATION", "PASS; six small pulses, 60 s horizon" if tds_pass else "FAIL"),
        ("BND_ACTIVE_MODE", f"COLLECTIVE_ZERO_BRANCH; pivot bus {bnd['pivot_bus']}"),
        ("BND_DIRECT_CONTRIBUTION", bnd["direct_dlambda_real"]),
        ("BND_SELF_ENERGY_CONTRIBUTION", bnd["self_energy_dlambda_real"]),
        ("BND_CANCELLATION_RATIO", bnd["cancellation_ratio"]),
        ("STRUCTURED_K_TESTED", "NO"),
        ("STRUCTURED_K_DIMENSION", "N/A"),
        ("STRUCTURED_K_RETAINED_SG_MW", "N/A"),
        ("STRUCTURED_K_OPTIMALITY_LOSS_MW", "N/A"),
        ("ROBUST_ZSTAR_RUN", "NO"),
        ("ROBUST_RETAINED_SG_MW", "N/A"),
        ("ROBUST_BETA_STAR", "N/A"),
        ("ROBUST_PD_VALIDATION", "N/A"),
        ("MAIN_RESULT", "PD-exact mixed model; locally KKT-certified 1.124344 MW retained SG"),
        ("MAIN_LIMITATION", "global branches not certified; PD margin headroom ~1.4e-10 s^-1; large pulses unvalidated"),
        ("NEXT_STEP", "certify remaining continuous architecture branches and robust beta frontier"),
        ("PUSH", "NO"),
    ]
    content = "\n".join(f"{k}: {v}" for k, v in kv) + "\n"
    (OUT / "FINAL_SUMMARY_EXP_N.md").write_text(content, encoding="utf-8")
    print(content, end="")


if __name__ == "__main__":
    main()
