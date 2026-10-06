"""Compare freshly re-executed IEEE-39 evidence with the frozen prior campaign."""

from pathlib import Path
import tomllib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "ieee39_exact_action_space_codesign_20261002"


def toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))


now = toml(HERE / "Q0_RESULT.toml")
old = toml(PRIOR / "Q0_RESULT.toml")
events = pd.read_csv(HERE / "Q0_EVENT_METRICS.csv")
old_events = pd.read_csv(PRIOR / "Q0_EVENT_METRICS.csv")
events = events.merge(old_events, on="event", suffixes=("_new", "_prior"), validate="one_to_one")
for col in ("F_peak_Hz", "RoCoF_peak_Hz_s", "Vmin_pu", "Vmax_pu", "min_SG_actuator_fraction_slack"):
    events[col + "_absolute_difference"] = (events[col + "_new"] - events[col + "_prior"]).abs()
events.to_csv(HERE / "M0_EVENT_REPRODUCTION.csv", index=False)
max_event_difference = max(events[c].max() for c in events if c.endswith("_absolute_difference"))

baseline = pd.DataFrame([{
    "case_id": "uniform_875_tau0",
    "GFL_percent": now["GFL_percent"],
    "GFL_MW": now["GFL_MW"],
    "retained_SG_MW": now["retained_SG_MW"],
    "equilibrium_residual_inf": now["equilibrium_residual_inf"],
    "physical_eigenvalues": now["physical_eigenvalues"],
    "critical_real_s_inv": now["critical_real_part_s_inv"],
    "GFL_percent_difference_vs_prior": abs(now["GFL_percent"] - old["GFL_percent"]),
    "retained_MW_difference_vs_prior": abs(now["retained_SG_MW"] - old["retained_SG_MW"]),
    "critical_alpha_difference_vs_prior": abs(now["critical_real_part_s_inv"] - old["critical_real_part_s_inv"]),
    "all_five_events_complete": bool(events.complete_new.all()),
    "all_five_events_pass": bool(events.pass_new.all()),
    "max_event_metric_difference_vs_prior": max_event_difference,
}])
baseline.to_csv(HERE / "M0_BASELINE_REPRODUCTION.csv", index=False)

q12 = toml(HERE / "Q1_Q2_STATUS.toml")
q12_old = toml(PRIOR / "Q1_Q2_STATUS.toml")
q4 = toml(HERE / "Q4_STATUS.toml")
q4_old = toml(PRIOR / "Q4_STATUS.toml")
q6 = toml(HERE / "Q6_STATUS.toml")
q6_old = toml(PRIOR / "Q6_STATUS.toml")
pi_new = pd.read_csv(HERE / "TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv")
pi_old = pd.read_csv(PRIOR / "TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv")
keys = ["bus", "tau_ms", "target_frequency_Hz"]
pi = pi_new.merge(pi_old, on=keys, suffixes=("_new", "_prior"), validate="one_to_one")
gain_difference = max(
    (pi[c + "_new"] - pi[c + "_prior"]).abs().max()
    for c in ("Kp_closed_form", "Ki_closed_form")
)
action = pd.DataFrame([{
    "case_id": "uniform_875_fixed_architecture",
    "PLL_delayed_channels": q12["delayed_channels"],
    "Q1_outer_product_residual": q12["max_Q1_relative_residual"],
    "Q2_operator_residual": q12["max_Q2_operator_relative_residual"],
    "Q2_operator_residual_difference_vs_prior": abs(q12["max_Q2_operator_relative_residual"] - q12_old["max_Q2_operator_relative_residual"]),
    "PI_cases": q4["tests"],
    "PI_full_characteristic_residual": q4["max_full_NEp_residual"],
    "PI_gain_max_absolute_difference_vs_prior": gain_difference,
    "Q6_action_dimension_upper_bound": q6["action_rank_upper_bound"],
    "Q6_factorization_residual": q6["max_factorization_relative_residual"],
    "Q6_factorization_residual_difference_vs_prior": abs(q6["max_factorization_relative_residual"] - q6_old["max_factorization_relative_residual"]),
}])
action.to_csv(HERE / "M0_ACTION_SPACE_REPRODUCTION.csv", index=False)
pi.to_csv(HERE / "M0_PI_REPRODUCTION.csv", index=False)

schur = pd.read_csv(HERE / "M0_ZERO_FREQUENCY_SCHUR_REPRODUCTION.csv").iloc[0]
pass_gate = (
    now["equilibrium_residual_inf"] < 1e-8
    and now["physical_eigenvalues"] == old["physical_eigenvalues"]
    and baseline.iloc[0]["critical_alpha_difference_vs_prior"] < 1e-6
    and baseline.iloc[0]["all_five_events_complete"]
    and baseline.iloc[0]["all_five_events_pass"]
    and max_event_difference < 1e-4
    and q12["max_Q2_operator_relative_residual"] < 1e-10
    and q4["max_full_NEp_residual"] < 1e-10
    and gain_difference < 1e-6
    and q6["max_factorization_relative_residual"] < 1e-10
    and schur["status"] == "BLOCKED_SINGULAR_ZERO_FREQUENCY_SCHUR"
)
lines = [
    "# M0 reproduction status",
    "",
    f"**M0 gate: {'PASS' if pass_gate else 'BLOCKED_BASELINE_REPRODUCTION'}.**",
    "",
    f"Fresh physical equilibrium residual: `{now['equilibrium_residual_inf']:.3e}`; physical poles: `{now['physical_eigenvalues']}`; critical real part: `{now['critical_real_part_s_inv']:.9f} s^-1`.",
    f"Five complete frozen events pass: `{bool(events.pass_new.all())}`. Largest difference from prior event metrics: `{max_event_difference:.3e}`. The copied validator re-executes the same frozen model; this is run parity, not a separate implementation of every physical component.",
    f"Low-rank operator reconstruction: `{q12['max_Q2_operator_relative_residual']:.3e}`; full descriptor action reconstruction: `{q6['max_factorization_relative_residual']:.3e}`. Conditional PI target cases: `{q4['tests']}` with characteristic residual `{q4['max_full_NEp_residual']:.3e}` and gain difference from prior `{gain_difference:.3e}`.",
    f"The separate gauge-deflated zero-frequency Schur recheck found hidden-block condition `{schur['hidden_condition']:.3e}` and status `{schur['status']}`. The invalid zero-frequency graph-damping construction remains excluded.",
    "",
    "The prior 90.047% joint design is not a five-event feasible baseline. M0 accepts only the re-executed 87.5% witness; no maximum is inferred.",
]
(HERE / "REPRODUCTION_STATUS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("M0_GATE", "PASS" if pass_gate else "BLOCKED_BASELINE_REPRODUCTION")
if not pass_gate:
    raise SystemExit(2)
