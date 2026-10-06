"""Print and save the evidence-qualified ExpK terminal summary."""
from pathlib import Path
import csv
import hashlib
import math
import tomllib

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_K"
TABLES = OUT / "tables"


def read_csv(name):
    path = TABLES / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def frozen(name):
    path = OUT / name
    if not path.exists():
        return None
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == (OUT / (name + ".sha256")).read_text().strip()
    return tomllib.loads(path.read_text(encoding="utf-8"))


nom = frozen("Z_K_NOMINAL_FINAL.toml")
rob = frozen("Z_K_ROBUST_beta_1p6991206999182038em6.toml")
catalog = read_csv("TABLE_K01_support_catalog.csv")
best = read_csv("TABLE_K08_support_best_candidates.csv")
front = read_csv("TABLE_K09_robust_frontier.csv")
trans = read_csv("TABLE_K14_transient_capacity.csv")
comparison = read_csv("TABLE_K11_blinded_comparison.csv")
bnd = read_csv("TABLE_K10_BND_explanation.csv")
graph = read_csv("TABLE_K10_graph_diagnostics.csv")
pd = read_csv("TABLE_K12_powerdynamics_validation.csv")
tds = read_csv("TABLE_K13_tds_validation.csv")


def value(obj, key, default="NOT_EVALUATED"):
    return obj.get(key, default) if obj else default


def pd_status(label):
    rows = [r for r in pd if r.get("candidate") == label]
    return rows[0]["status"] + " alpha_PD=" + rows[0]["alpha_PD"] if rows else "NOT_RUN"


def tds_status(label):
    rows = [r for r in tds if r.get("candidate") == label]
    if not rows:
        return "NOT_RUN"
    scaling = all(r["status"] == "EVALUATED" and r.get("finite", "").lower() == "true" and
                         float(r.get("relative_frequency_scaling_error", "nan")) <= .1 and
                         float(r.get("relative_voltage_scaling_error", "nan")) <= .1
                         for r in rows)
    if pd_status(label).startswith("FAIL"):
        return "SCALING_ONLY_PD_FAILED" if scaling else "FAIL_PD_AND_TDS"
    return "PASS" if scaling else "FAIL"


nomtrans = sorted((r for r in trans if r["candidate"] == "nominal"),
                  key=lambda r: float(r["deltaP_max_total_MW"]))
worst = nomtrans[0] if nomtrans else {}
erow = next((r for r in comparison if r["case"] == "ExpE_provisional"), {})
grow = next((r for r in comparison if r["case"] == "ExpG_final"), {})
psi = [float(r["Psi_MW"]) for r in read_csv("TABLE_K06_dual_survival_values.csv")
       if r["status"] == "DIAGNOSTIC_MULTIPLIER_NOT_KKT_CERTIFIED"]
kkt_diag = math.sqrt(sum(x*x for x in psi)) if psi else float("nan")
pdnom = pd_status("nominal")
pdrob = pd_status("robust_ExpG_beta")
tdsnom = tds_status("nominal")
tdsrob = tds_status("robust_ExpG_beta")
status = "FAIL" if pdnom.startswith("FAIL") or pdrob.startswith("FAIL") else "PARTIAL"
lines = [
    f"EXP_K_STATUS: {status}",
    f"TOTAL_SG_MW: {sum(nom['P_initial_MW']) if nom else 'NOT_EVALUATED'}",
    f"SUPPORTS_EVALUATED: {len(catalog)}",
    "SUPPORTS_CERTIFIED_INFEASIBLE: 0",
    f"SUPPORTS_WITH_FEASIBLE_BRANCHES: {sum(r['status'] != 'NO_FEASIBLE_UNIFORM_SEED_NOT_PROVEN_INFEASIBLE' for r in best)} sampled seeds",
    "",
    f"NOMINAL_BEST_SUPPORT: {value(nom, 'support_buses')}",
    f"NOMINAL_RHO: {value(nom, 'rho')}",
    f"NOMINAL_KP: {value(nom, 'Kp')}",
    f"NOMINAL_KI: {value(nom, 'Ki')}",
    f"NOMINAL_RETAINED_SG_MW: {value(nom, 'retained_SG_MW')}",
    f"NOMINAL_GFL_MW: {value(nom, 'GFL_MW')}",
    f"NOMINAL_GFL_FRACTION: {nom['GFL_MW']/sum(nom['P_initial_MW']) if nom else 'NOT_EVALUATED'}",
    f"NOMINAL_ALPHA: {value(nom, 'spectral_abscissa_s_inv')}",
    f"NOMINAL_ACTIVE_POLES: {list(zip(nom['active_poles_real'],nom['active_poles_imag'])) if nom else 'NOT_EVALUATED'}",
    f"NOMINAL_KKT_RESIDUAL: NOT_CERTIFIED; interior diagnostic={kkt_diag} MW",
    "NOMINAL_SOSC: NOT_EVALUATED",
    "NOMINAL_GLOBAL_CERTIFIED: NO",
    "",
    f"ROBUST_FRONTIER: {len(front)} normalized-beta points; BEST_SAMPLED_ONLY",
    "ROBUST_BETA_SELECTED: 1.6991206999182038e-6 (ExpG normalized comparison)",
    f"ROBUST_BEST_SUPPORT: {value(rob, 'support_buses')}",
    f"ROBUST_RHO: {value(rob, 'rho')}",
    f"ROBUST_KP: {value(rob, 'Kp')}",
    f"ROBUST_KI: {value(rob, 'Ki')}",
    f"ROBUST_RETAINED_SG_MW: {value(rob, 'retained_SG_MW')}",
    f"ROBUST_GFL_MW: {value(rob, 'GFL_MW')}",
    f"ROBUST_GFL_FRACTION: {rob['GFL_MW']/sum(rob['P_initial_MW']) if rob else 'NOT_EVALUATED'}",
    f"ROBUST_ALPHA: {value(rob, 'spectral_abscissa_s_inv')}",
    f"ROBUST_BETA_STAR: {value(rob, 'beta_star')}",
    f"ROBUST_SMALL_GAIN_PEAK: {value(rob, 'small_gain_peak')}",
    f"ROBUST_PEAK_FREQUENCY: {value(rob, 'peak_frequency_rad_s')}",
    "ROBUST_GLOBAL_CERTIFIED: NO",
    "",
    f"WORST_DISTURBANCE_BUS: {value(worst, 'bus')}",
    f"ROCOF_UNIT_GAIN: {value(worst, 'peak_RoCoF_Hz_s_per_MW')}",
    f"FREQUENCY_UNIT_GAIN: {value(worst, 'peak_frequency_Hz_per_MW')}",
    f"DELTA_P_MAX_ROCOF: {value(worst, 'deltaP_max_RoCoF_MW')}",
    f"DELTA_P_MAX_FREQUENCY: {value(worst, 'deltaP_max_frequency_MW')}",
    f"DELTA_P_MAX_TOTAL: {value(worst, 'deltaP_max_total_MW')}",
    "",
    f"EXP_E_NOMINAL_COMPARISON: E retained={value(erow,'retained_SG_MW')} MW; E alpha={value(erow,'exact_ExpK_alpha_s_inv')}; E BEATS K",
    f"EXP_G_EQUAL_BETA_COMPARISON: G GFL={value(grow,'GFL_MW')} MW; K robust GFL={value(rob,'GFL_MW')} MW analytically; PD result={pdrob}",
    f"BND_DOMINANT_ACTIVE_MODE: {value(bnd[0] if bnd else None,'pole_real')} + j{value(bnd[0] if bnd else None,'pole_imag')}",
    f"BND_SELF_ENERGY_SHARE: {value(bnd[0] if bnd else None,'self_energy_share')}",
    f"GRAPH_DIAGNOSTIC: L_G/L_B commutator ratio={value(graph[0] if graph else None,'LG_LB_commutator_ratio')}; no modal law claimed",
    "",
    f"POWERDYNAMICS_NOMINAL_VALIDATION: {pdnom}",
    f"POWERDYNAMICS_ROBUST_VALIDATION: {pdrob}",
    f"TDS_NOMINAL_VALIDATION: {tdsnom}",
    f"TDS_ROBUST_VALIDATION: {tdsrob}",
    "GENERIC_SOLVER_FALSIFICATION: NOT_RUN",
    "BETTER_POINT_FOUND_AFTER_FREEZE: YES (ExpE archived analytical point)",
    "MAIN_THEORETICAL_RESULT: Support transitions change state dimension; complete-spectrum and self-energy cancellation govern local stability.",
    "MAIN_POWER_SYSTEM_RESULT: Sampled near-full-replacement points fail independent PowerDynamics validation and have low linear disturbance capacity.",
    "MAIN_LIMITATION: No branch-complete KKT/global certificate; ExpE is better; independent validation fails.",
    "PUSH: NO",
]
summary = "\n".join(lines) + "\n"
(OUT / "FINAL_SUMMARY_EXP_K.md").write_text(summary, encoding="utf-8")
print(summary, end="")
