"""Create canonical top-level tables and explicit BLOCKED records; never impute results."""
from pathlib import Path
import csv
import shutil
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE / "baseline_reproduction"

for name in (
    "TABLE_D00_BASELINE_REPRODUCTION.csv",
    "TABLE_D01_EXACT_DDE_POLES.csv",
    "TABLE_D02_DERIVATIVE_VALIDATION.csv",
    "TABLE_D04_GRAPH_DAMPING_OPERATORS.csv",
):
    source = BASE / name
    if source.exists():
        shutil.copy2(source, HERE / name)

# Recompute graph-signal descriptors from the already-frozen vectors and graph
# matrices. Do not alter or re-rank the delay assignments.
lp = pd.read_csv(HERE / "GRAPH_LP_LOSSLESS_PORT.csv").to_numpy(float)
mass = pd.read_csv(HERE / "GRAPH_MASS_REFERENCE.csv")["M"].to_numpy(float)
U = pd.read_csv(HERE / "GRAPH_GSP_BASIS.csv").to_numpy(float)
eig = pd.read_csv(HERE / "GRAPH_MODES.csv")["eigenvalue"].to_numpy(float)
lbar = np.diag(1/np.sqrt(mass)) @ lp @ np.diag(1/np.sqrt(mass))
base_ms = np.arange(10, dtype=float) * (50/9)
frozen = pd.read_csv(HERE / "FROZEN_DELAY_PATTERNS.csv")
delay_cases = [(f"uniform_{t}ms", np.full(10, float(t)), "uniform")
               for t in (0,2,5,10,20,30,40,50)]
delay_cases += [(str(row.delay_pattern_id), np.array([float(x) for x in row.tau_ms.split(";")]), "heterogeneous")
                for row in frozen.itertuples(index=False)]
gsp_rows = []
for pattern_id, tau_ms, family in delay_cases:
    tau = tau_ms / 1000
    comm = lbar @ np.diag(tau) - np.diag(tau) @ lbar
    tauhat = U.T @ tau
    comm_physical = lp @ np.diag(tau) - np.diag(tau) @ lp
    comm_sq = np.linalg.norm(comm_physical,"fro")**2
    edge_identity_sq = 2*sum(lp[i,j]**2*(tau[i]-tau[j])**2
                             for i in range(10) for j in range(i+1,10))
    gsp_rows.append({
        "delay_pattern_id": pattern_id, "family": family,
        "mean_tau_ms": tau_ms.mean(), "std_tau_ms": tau_ms.std(), "max_tau_ms": tau_ms.max(),
        "chi_tau": np.linalg.norm(comm, "fro")/max(np.linalg.norm(lbar,"fro"), np.finfo(float).eps),
        "graph_roughness": tau @ lp @ tau,
        "tau_hat_seconds": ";".join(f"{x:.17g}" for x in tauhat),
        "commutator_identity_relative_error": abs(comm_sq-edge_identity_sq)/max(comm_sq,np.finfo(float).eps),
        "low_frequency_score": np.sum(np.maximum(eig,0)*tauhat**2),
        "high_frequency_score": np.sum(np.maximum(eig,0)**2*tauhat**2),
        "vector_sha256": "uniform_grid" if family=="uniform" else str(frozen.loc[frozen.delay_pattern_id==pattern_id,"vector_sha256"].iloc[0]),
        "graph_metrics_status": "NUMERICALLY_VALIDATED",
        "operator_status": "BLOCKED_NO_VALID_SCHUR_SLOPE"})
pd.DataFrame(gsp_rows).to_csv(HERE / "TABLE_D05_GSP_DELAY_METRICS.csv", index=False)

def write_rows(filename, fields, rows):
    with (HERE / filename).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

write_rows("TABLE_D06_GSP_GAIN_REDUCTION.csv",
    ["q", "max_GFL_MW", "max_GFL_percent", "compute_seconds", "status", "reason"],
    [{"q": q, "status": "BLOCKED_BY_D1", "reason": "No complete DDE spectrum; gain-restricted co-design not run"}
     for q in [0,1,2,3,4,5,"full_nodal"]])

d07_fields = ["delay_pattern_id","mean_tau_ms","std_tau_ms","chi_tau","graph_roughness","eta_cross",
 "max_GFL_MW","max_GFL_percent","retained_SG_MW","rho_i_star","Kp_i_star","Ki_i_star",
 "critical_mode","critical_frequency","critical_alpha","active_constraint","RoCoF",
 "frequency_excursion","full_event_feasible","status","reason"]
write_rows("TABLE_D07_REPLACEMENT_FRONTIER.csv", d07_fields,
    [{"status":"BLOCKED_EXACT_DDE_SPECTRUM","reason":"No complete rightmost-root coverage; no delay-dependent optimization or feasible maximum is claimed"}])

write_rows("TABLE_D08_DELAY_SHADOW_PRICES.csv",
    ["delay_pattern_id","bus","tau_ms","shadow_price_MW_per_ms","finite_change_MW_per_ms","active_set_stable","status","reason"],
    [{"status":"BLOCKED_BY_D1","reason":"No validated delayed optimum/KKT active set"}])

headline = ["tau0_baseline","uniform_delay","graph_smooth","graph_rough","best_validated_replacement","adverse_delay_perturbation"]
write_rows("TABLE_D09_NONLINEAR_DDE_VALIDATION.csv",
    ["case_id","event_id","max_frequency_Hz","max_RoCoF_Hz_s","settling_s","actuator_margin","guards_pass","status","reason"],
    [{"case_id": c, "event_id": "all_frozen_events", "status": "BLOCKED_NONLINEAR_DDE",
      "reason": "No reliable method-of-steps DDE integrator was established; delayed event simulation not run"} for c in headline])

write_rows("TABLE_D10_OPTIMALITY_BOUNDS.csv",
    ["delay_pattern_id","R_F_MW","R_U_MW","gap_percent","remainder_bound","bound_method","optimality_status","reason"],
    [{"optimality_status":"BLOCKED_NUMERICAL_ONLY_NOT_AVAILABLE",
      "reason":"No validated delayed feasible lower point and no rigorous Taylor remainder/global upper bound"}])

print("canonical tables copied; D06-D10 explicitly record unavailable outcomes")
