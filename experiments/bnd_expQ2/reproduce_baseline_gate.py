"""Read-only ExpN/ExpP/ExpQ baseline regression gate for ExpQ2."""
from __future__ import annotations
import csv, hashlib, json, tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_Q2" / "BASELINE"
OUT.mkdir(parents=True, exist_ok=True)
model_sha = "e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a"
q0 = json.loads((ROOT / "reports/experiment_Q/Q0/Q0_RESULTS.json").read_text())
all_sg = json.loads((ROOT / "reports/experiment_Q/Q2_REFERENCE/Q2_ALL_SG_REFERENCE.json").read_text())
expn_path = ROOT / "reports/experiment_N/Z_N_NOMINAL_FINAL.toml"
expp_path = ROOT / "reports/experiment_P/P5/Z_P_NOMINAL_FINAL.toml"
expn = tomllib.loads(expn_path.read_text())
expp = tomllib.loads(expp_path.read_text())
row_checks = []

def check(name, value, expected, passed):
    row_checks.append({"check": name, "value": value, "expected": expected,
                       "pass": bool(passed)})

check("model_sha", q0["model_sha"], model_sha, q0["model_sha"] == model_sha)
for label, path, expected in (("ExpN candidate SHA", expn_path, q0["ExpN_candidate_sha256"]),
                              ("ExpP candidate SHA", expp_path, q0["ExpP_candidate_sha256"])):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    check(label, digest, expected, digest == expected)
check("ExpN support bus", expn["support_bus"], 38, expn["support_bus"] == 38)
check("ExpN retained MW", expn["retained_SG_MW"], 1.124344477653, abs(expn["retained_SG_MW"]-1.124344477653)<1e-8)
check("ExpN rho38", expn["rho"][8], 0.9986453680992127, abs(expn["rho"][8]-0.9986453680992127)<1e-12)
check("ExpN alpha", q0["alpha_analytic"], -0.050000001105, abs(q0["alpha_analytic"]+0.050000001105)<1e-8)
check("ExpN PD alpha", q0["alpha_PD"], -0.050000000139, abs(q0["alpha_PD"]+0.050000000139)<1e-8)
check("ExpN robust infeasible witness", q0["pointwise_beta_upper"], q0["beta_req"], q0["pointwise_beta_upper"] < q0["beta_req"])
check("P/Q trim", max(q0["max_P_error_pu"], q0["max_Q_error_pu"]), 1e-12,
      max(q0["max_P_error_pu"], q0["max_Q_error_pu"]) < 1e-12)
check("ExpN 100MW step frequency failure", q0["step_peak_frequency_Hz"], 39.3334,
      q0["step_peak_frequency_Hz"] > 0.5)
check("ExpN 100MW step RoCoF failure", q0["step_peak_RoCoF_Hz_s"], 1.07221,
      q0["step_peak_RoCoF_Hz_s"] > 0.5)
check("all-SG event feasible", all_sg["reference_feasible"], True,
      all_sg["reference_feasible"] and all_sg["PD_F_peak"] < 0.5 and all_sg["PD_R_peak"] < 0.5)
check("all-SG full spectrum margin", all_sg["alpha"], -0.05, all_sg["alpha"] <= -0.05)
status = "PASS_BASELINE" if all(r["pass"] for r in row_checks) else "FAIL_BASELINE"
with (OUT / "TABLE_Q2_BASELINE_REGRESSION.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["check", "value", "expected", "pass"])
    writer.writeheader(); writer.writerows(row_checks)
summary = {
    "status": status, "model_sha": model_sha, "checks": len(row_checks),
    "all_pass": all(r["pass"] for r in row_checks),
    "ExpN_alpha_analytic": q0["alpha_analytic"], "ExpN_alpha_PD": q0["alpha_PD"],
    "ExpN_retained_SG_MW": q0["retained_SG_MW"],
    "ExpN_rho38": expn["rho"][8], "ExpN_nominal_step_peak_Hz": q0["step_peak_frequency_Hz"],
    "ExpN_nominal_step_RoCoF_Hz_s": q0["step_peak_RoCoF_Hz_s"],
    "ExpN_beta_pointwise_upper": q0["pointwise_beta_upper"],
    "max_PQ_error_pu": max(q0["max_P_error_pu"], q0["max_Q_error_pu"]),
    "allSG_alpha": all_sg["alpha"], "allSG_PD_frequency_peak_Hz": all_sg["PD_F_peak"],
    "allSG_PD_RoCoF_Hz_s": all_sg["PD_R_peak"],
}
(OUT / "BASELINE_RESULTS.json").write_text(json.dumps(summary, indent=2) + "\n")
(OUT / "STAGE_SUMMARY.md").write_text(f"""# ExpQ2 baseline regression gate

**Status: {status}.** This gate reuses the just-reproduced ExpQ Q0 PowerDynamics build/TDS results and the frozen all-SG PowerDynamics TDS reference, then independently checks their expected values and the immutable ExpN/ExpP candidate hashes. It makes no design-model or candidate changes.

- Passed checks: {sum(x['pass'] for x in row_checks)}/{len(row_checks)}.
- ExpN: support 38, retained {q0['retained_SG_MW']:.12f} MW, rho38 {expn['rho'][8]:.15f}, alpha analytic/PD {q0['alpha_analytic']:.12f}/{q0['alpha_PD']:.12f} s⁻¹.
- P/Q max error {summary['max_PQ_error_pu']:.3e} pu; ExpN 100 MW bus-16 peak/RoCoF {q0['step_peak_frequency_Hz']:.6f} Hz / {q0['step_peak_RoCoF_Hz_s']:.6f} Hz/s; pointwise beta witness {q0['pointwise_beta_upper']:.4e} < {q0['beta_req']:.4e}.
- All-SG PD TDS peak/RoCoF {all_sg['PD_F_peak']:.6f} Hz / {all_sg['PD_R_peak']:.6f} Hz/s, alpha {all_sg['alpha']:.9f} s⁻¹.
- Evidence: `TABLE_Q2_BASELINE_REGRESSION.csv`; no optimizer has run.
""", encoding="utf-8")
print(f"{status}: {sum(x['pass'] for x in row_checks)}/{len(row_checks)} baseline checks")
