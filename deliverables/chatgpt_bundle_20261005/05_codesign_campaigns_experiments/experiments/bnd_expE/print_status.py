"""Print the preregistered ExpE terminal summary without implying PASS."""
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
r = json.loads((root / "reports" / "experiment_E" / "RESULTS_EXP_E.json").read_text(encoding="utf-8"))
values = [
    ("EXP_E_STATUS", "PARTIAL"),
    ("DESIGN_USES_POWERDYNAMICS", "NO"),
    ("DESIGN_USES_GENERIC_OPTIMIZER", "NO"),
    ("REPLACEABLE_GENERATOR_COUNT", r["replaceable_generator_count"]),
    ("REPLACEABLE_BUSES", ",".join(map(str, r["replaceable_buses"]))),
    ("DESIGN_VARIABLE_COUNT", r["design_variable_count"]),
    ("PLL_KP_KI_INDEPENDENT", "YES"),
    ("TOTAL_INITIAL_SG_MW", r["actual_initial_SG_MW"]),
    ("PRIMARY_OPTIMALITY_STATUS", "BEST_ANALYTIC_FEASIBLE"),
    ("TOTAL_MAX_GFL_MW", r["provisional_GFL_MW"]),
    ("TOTAL_REMAINING_SG_MW", r["provisional_retained_SG_MW"]),
    ("REPLACEMENT_PERCENT", 100*r["provisional_fraction"]),
    ("STATIC_ENDPOINT_FEASIBLE", "NO"),
    ("CONTINUOUS_PATH_STATUS", "INCONCLUSIVE"),
    ("MIN_EFFORT_GAINS", "NOT_COMPUTED"),
    ("MAX_MARGIN_GAINS", "NOT_COMPUTED"),
    ("MINIMUM_FINAL_SPECTRAL_MARGIN", -r["provisional_spectral_abscissa_per_s"]),
    ("MINIMUM_FINAL_DAMPING_RATIO", "NOT_COMPUTED"),
    ("ACTIVE_LIMITING_MODES", "real pole near -0.05 1/s; angle gauge excluded"),
    ("PLL_AUTHORITY_MIN_SINGULAR_VALUE", "NOT_CERTIFIED_AT_FINAL_POINT"),
    ("COLLECTIVE_NONCOMPOSABILITY_FOUND", "YES"),
    ("STABILITY_HOLE_FOUND", "NO"),
    ("SELF_ENERGY_MECHANISM", "INCONCLUSIVE"),
    ("ANALYTIC_CANDIDATE_SHA256", "NOT_FROZEN"),
    ("POWERDYNAMICS_VALIDATION", "NOT_RUN"),
    ("ANALYTIC_VS_PD_MAX_POLE_ERROR", "NOT_AVAILABLE"),
    ("FINAL_TABLE", "NOT_AVAILABLE"),
    ("MAIN_ANALYTIC_RESULT", "A 20-port closure yields a strictly feasible provisional 5401.571878 MW allocation."),
    ("MAIN_POWER_SYSTEM_RESULT", "Initialized SG dispatch is 5402.761090 MW; all-GFL has a near-double-zero frequency mode."),
    ("MAIN_BND_RESULT", "The all-SG and four single-bus spectra reproduce frozen ExpC within 1.76e-11."),
    ("MAIN_LIMITATION", "Global/branch optimality, final freeze, independent detailed-model and TDS validation remain open."),
    ("PUSH", "NO"),
]
for label, value in values:
    print(f"{label}:\n{value}\n")
