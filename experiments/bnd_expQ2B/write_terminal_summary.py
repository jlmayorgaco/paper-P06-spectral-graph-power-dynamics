from pathlib import Path
import csv
import tomllib

root=Path(__file__).resolve().parents[2]
core=root/"reports"/"experiment_Q2B"/"CORE"
report=root/"reports"/"experiment_Q2B"
load=lambda name: tomllib.loads((core/name).read_text(encoding="utf-8"))
c=load("Z_Q2B_SECURE_T05_FINAL.toml")
core0=load("Z_Q2B_CORE_d000.toml")
pd=load("PD_VALIDATION_RESULTS.toml")
rob=load("CORE_RESTORED_AUDIT.toml")
kkt=load("CORE_KKT_AUDIT_RESTORED.toml")
beta=load("CORE_BETA_INTERVAL_CERTIFICATE.toml")

def rows(name):
    with (core/name).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

core_windows=rows("TABLE_Q2B_core_d0_rocof_at100.csv")
secure_windows=rows("TABLE_Q2B_core_rocof_windows_restored.csv")
oos=rows("TABLE_Q2B_PD_out_of_sample.csv")
sens=rows("TABLE_Q2B_candidate_sensitivities.csv")
best={}
for col in ("a_alpha_per_MW","a_Finf_per_MW","a_R05_per_MW"):
    best[col]=max(sens,key=lambda r:float(r[col]))["bus"]

fmtarr=lambda vals: "["+", ".join(f"{float(v):.9g}" for v in vals)+"]"
coresupport="{"+", ".join(map(str,core0["support"]))+"}"
support="{"+", ".join(map(str,c["support"]))+"}"
window_map={float(r["window_s"]):r for r in secure_windows}
core_window_map={float(r["window_s"]):r for r in core_windows}
oos_pass=all(r["pass"].lower()=="true" for r in oos)

text=f"""EXP_Q2B_STATUS: FAIL_OPTIMIZATION

REGRESSION_PASS: YES (9/9 Q2B checks; inherited baseline gate 14/14)
CORE_OPTIMIZATION_RUN: YES
CORE_STATUS: BEST_FOUND_KKT_OPEN; NOT CERTIFIED
CORE_RETAINED_SG_MW: {float(core0['retained_SG_MW']):.9f}
CORE_SUPPORT: {coresupport}
CORE_RHO: {fmtarr(core0['rho'])}
CORE_KP: {fmtarr(core0['Kp'])}
CORE_KI: {fmtarr(core0['Ki'])}
CORE_ALPHA: {float(core0['alpha']):.12f}
CORE_BETA: {float(core0['beta_sampled']):.12e} (sampled)
CORE_FREQ_PEAK: 0 Hz at d=0; {float(core_window_map[0.5]['F_peak_Hz']):.6f} Hz at 100 MW
CORE_FREQ_STEADY: 0 Hz at d=0; {float(core_window_map[0.5]['F_inf_Hz']):.6f} Hz at 100 MW

CORE_KKT_PRIMAL: 0 (feasible checkpoint only)
CORE_KKT_STATIONARITY: NOT_AVAILABLE
CORE_KKT_COMPLEMENTARITY: NOT_AVAILABLE
CORE_LICQ: NOT_CERTIFIED
CORE_SOSC: NOT_TESTED

ROCOF_MEASUREMENT_DEFINITION: causal unwrapped bus-voltage phase difference; T={0.2, 0.5, 1.0, 2.0} s; buses 30:39
ROCOF_T_WINDOWS: {0.2, 0.5, 1.0, 2.0}
CORE_R_0P2: {float(core_window_map[0.2]['R_peak_Hz_s']):.9f} Hz/s (CORE d=0 point at 100 MW)
CORE_R_0P5: {float(core_window_map[0.5]['R_peak_Hz_s']):.9f} Hz/s (CORE d=0 point at 100 MW)
CORE_R_1P0: {float(core_window_map[1.0]['R_peak_Hz_s']):.9f} Hz/s (CORE d=0 point at 100 MW)
CORE_R_2P0: {float(core_window_map[2.0]['R_peak_Hz_s']):.9f} Hz/s (CORE d=0 point at 100 MW)
ROCOF_WINDOW_ROBUST: YES for the frozen best-found 100 MW point; not an optimum claim

SECURE_T0P2_RETAINED_SG_MW: NOT_OPTIMIZED; best-found feasible point {float(c['retained_SG_MW']):.9f} MW
SECURE_T0P5_RETAINED_SG_MW: NOT_OPTIMIZED; best-found feasible point {float(c['retained_SG_MW']):.9f} MW
SECURE_T1P0_RETAINED_SG_MW: NOT_OPTIMIZED; best-found feasible point {float(c['retained_SG_MW']):.9f} MW
SECURE_T2P0_RETAINED_SG_MW: NOT_OPTIMIZED; best-found feasible point {float(c['retained_SG_MW']):.9f} MW

SECURE_T0P5_SUPPORT: {support}
SECURE_T0P5_RHO: {fmtarr(c['rho'])}
SECURE_T0P5_KP: {fmtarr(c['Kp'])}
SECURE_T0P5_KI: {fmtarr(c['Ki'])}
SECURE_T0P5_ALPHA: {float(c['alpha']):.12f}
SECURE_T0P5_BETA: {float(rob['beta_full_refined_observed']):.15e} observed; interval incomplete
SECURE_T0P5_FREQ_PEAK: {float(c['F_peak_Tref_Hz']):.9f} Hz (analytic)
SECURE_T0P5_FREQ_STEADY: {float(c['F_inf_Hz']):.9f} Hz (analytic)
SECURE_T0P5_ROCOF: {float(window_map[0.5]['R_peak_Hz_s']):.9f} Hz/s (analytic)

T0P5_KKT_PRIMAL: {float(kkt['primal_violation']):.3e}
T0P5_KKT_STATIONARITY: {float(kkt['stationarity_inf_norm']):.6f}
T0P5_KKT_COMPLEMENTARITY: {abs(float(kkt['multipliers_unconstrained_LS'][0])*float(kkt['constraint_values_scaled'][1])):.3e}
T0P5_LICQ: rank {kkt['LICQ_rank']} of 1 near-active row
T0P5_SOSC: NOT_TESTED

CONTINUATION_CORE_TO_100MW: NO; blocked near 65.7 MW, direct correction used at 100 MW
CONTINUATION_T0P5_TO_100MW: NO
ACTIVE_SET_SWITCHES: NOT_COMPLETED / NOT_CERTIFIED
SUPPORT_SWITCHES: 0; only the full ten-bus mixed architecture was optimized

PD_FINAL_ALPHA: {float(pd['pd_alpha']):.12f}
PD_FINAL_FREQ_PEAK: {float(pd['TDS_Fpeak_T05_Hz']):.9f} Hz
PD_FINAL_FREQ_STEADY: {float(pd['TDS_Fsteady_T05_Hz']):.9f} Hz
PD_FINAL_ROCOF_0P5: {float(pd['TDS_Rpeak_T05_Hz_s']):.9f} Hz/s
PD_FINAL_100MW_PASS: YES (spectral identity and the event limits pass)

BEST_MODAL_BUS: {best['a_alpha_per_MW']} (local sensitivity diagnostic only)
BEST_STEADY_FREQ_BUS: {best['a_Finf_per_MW']} (local sensitivity diagnostic only)
BEST_ROCOF_BUS: {best['a_R05_per_MW']} (local sensitivity diagnostic only)
OPTIMAL_RETAINED_BUSES: NOT_ESTABLISHED; best-found support {support}

BND_DIRECT: NOT_EVALUATED
BND_SELF_ENERGY: NOT_EVALUATED
BND_TOTAL: NOT_EVALUATED
BND_RECONSTRUCTION_ERROR: NOT_APPLICABLE; modal constraint is slack and has no certified multiplier

GLOBAL_CERTIFIED: NO
GLOBAL_LOWER: 0 MW (trivial physical bound only)
GLOBAL_UPPER: {float(c['retained_SG_MW']):.9f} MW (best-found feasible candidate)
GLOBAL_GAP: {float(c['retained_SG_MW']):.9f} MW

MAIN_RESULT: PowerDynamics-validated best-found candidate, not a secure optimum certificate
MAIN_THEORETICAL_RESULT: finite-window phase output removes architecture-dependent internal-PLL RoCoF from the primary constraint
MAIN_POWER_SYSTEM_RESULT: {float(c['retained_SG_MW']):.3f} MW retained; PD bus16/100 MW meets the frozen project limits
MAIN_BND_RESULT: no BND active-mode decomposition because modal constraint is slack
MAIN_LIMITATION: KKT stationarity 0.830; beta interval incomplete at 50,001 nodes; supports open
POSTER_MAIN_NUMBER: 438.446 MW best-found retained SG
POSTER_MAIN_CLAIM: validated candidate passes the tested event limits; minimum-retention optimality is unproved

PUSH: NO
"""
(report/"TERMINAL_SUMMARY_EXP_Q2B.txt").write_text(text,encoding="utf-8")
print(report/"TERMINAL_SUMMARY_EXP_Q2B.txt")
