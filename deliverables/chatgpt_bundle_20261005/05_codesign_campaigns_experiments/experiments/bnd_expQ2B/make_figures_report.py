from __future__ import annotations

import csv
import hashlib
import math
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports" / "experiment_Q2B"
CORE = REPORT / "CORE"
FIG = REPORT / "FIGURES"
FIG.mkdir(parents=True, exist_ok=True)

with (CORE / "Z_Q2B_SECURE_T05_FINAL.toml").open("rb") as f:
    cand = tomllib.load(f)
with (CORE / "CORE_RESTORED_AUDIT.toml").open("rb") as f:
    robust_audit = tomllib.load(f)
with (CORE / "CORE_KKT_AUDIT_RESTORED.toml").open("rb") as f:
    kkt = tomllib.load(f)
with (CORE / "CORE_BETA_INTERVAL_CERTIFICATE.toml").open("rb") as f:
    beta_interval = tomllib.load(f)
with (CORE / "PD_VALIDATION_RESULTS.toml").open("rb") as f:
    pd = tomllib.load(f)

candidate_bus = pandas.read_csv(CORE / "TABLE_Q2B_candidate_by_bus.csv")
windows_an = pandas.read_csv(CORE / "TABLE_Q2B_core_rocof_windows_restored.csv")
core_d0_windows = pandas.read_csv(CORE / "TABLE_Q2B_core_d0_rocof_at100.csv")
windows_pd = pandas.read_csv(CORE / "TABLE_Q2B_PD_TDS_metrics.csv")
continuation = pandas.read_csv(CORE / "TABLE_Q2B_core_continuation.csv")
oos = pandas.read_csv(CORE / "TABLE_Q2B_PD_out_of_sample.csv")
sensitivity = pandas.read_csv(CORE / "TABLE_Q2B_candidate_sensitivities.csv")
analytic_trace = pandas.read_csv(CORE / "TABLE_Q2B_analytic_TDS_timeseries.csv")
pd_trace = pandas.read_csv(CORE / "TABLE_Q2B_PD_TDS_timeseries.csv")

Ptotal = float(candidate_bus.Pgen0_MW.sum())
J = float(cand["retained_SG_MW"])
GFL = Ptotal - J
fraction = GFL / Ptotal
beta_req = float(cand["beta_req_normalized"])
beta_obs = float(robust_audit["beta_full_refined_observed"])
alpha = float(cand["alpha"])
finf = float(cand["F_inf_Hz"])
fpeak = float(cand["F_peak_Tref_Hz"])
rpeak = float(pd["TDS_Rpeak_T05_Hz_s"])

# Consolidated candidate record. The full-frequency result is observed, while
# the interval routine remained incomplete; do not promote its partial bound.
constraint_rows = [
    ("modal alpha", alpha, -0.05, (-0.05-alpha)/0.05, "OBSERVED_PASS"),
    ("robust beta-star", beta_obs, beta_req, (beta_obs-beta_req)/beta_req,
     "OBSERVED_MARGIN; INTERVAL_INCOMPLETE"),
    ("steady grid frequency", finf, 0.5, (0.5-finf)/0.5, "ANALYTIC_PASS"),
    ("peak grid frequency T=0.5", fpeak, 0.5, (0.5-fpeak)/0.5, "ANALYTIC_PASS"),
    ("RoCoF T=0.5", float(windows_an.loc[windows_an.window_s == 0.5, "R_peak_Hz_s"].iloc[0]),
     0.5, (0.5-float(windows_an.loc[windows_an.window_s == 0.5, "R_peak_Hz_s"].iloc[0]))/0.5,
     "ANALYTIC_PASS"),
]
pandas.DataFrame(constraint_rows, columns=["constraint", "observed_value", "limit",
    "normalized_slack", "status"]).to_csv(CORE / "TABLE_Q2B_frozen_candidate_constraints.csv", index=False)

# Record final stage rows while keeping the original attempted continuation log.
cont_final = continuation[continuation.disturbance_MW < 100.0].copy()
final_row = pandas.DataFrame([{
    "disturbance_MW": 100.0, "status": "BEST_FOUND_NOT_LOCAL_KKT",
    "cost": J, "alpha": alpha, "beta": beta_obs, "Finf": finf,
    "Fpeak": fpeak, "max_g_sampled": max(kkt["constraint_values_scaled"]),
    "iterations": 6, "accepted_steps": 1,
    "KKT_stationarity": float(kkt["stationarity_inf_norm"]),
    "KKT_primal": float(kkt["primal_violation"]),
}])
cont_final["max_g_sampled"] = cont_final["max_g"]
cont_final = cont_final.drop(columns=["max_g"])
pandas.concat([cont_final, final_row], ignore_index=True).to_csv(
    CORE / "TABLE_Q2B_continuation_and_best_candidate.csv", index=False)
final_metrics = [
    ("status", "FAIL_OPTIMIZATION"), ("model_sha", cand["model_sha"]),
    ("candidate_sha256", hashlib.sha256((CORE / "Z_Q2B_SECURE_T05_FINAL.toml").read_bytes()).hexdigest()),
    ("total_original_SG_MW", Ptotal), ("best_found_retained_SG_MW", J),
    ("best_found_GFL_MW", GFL), ("GFL_fraction", fraction),
    ("alpha_analytic_s_inv", alpha), ("beta_full_refined_observed", beta_obs),
    ("beta_interval_status", beta_interval["status"]),
    ("F_inf_analytic_Hz", finf), ("F_peak_T05_analytic_Hz", fpeak),
    ("R_peak_T05_analytic_Hz_s", float(windows_an.loc[windows_an.window_s==.5,"R_peak_Hz_s"].iloc[0])),
    ("KKT_stationarity", float(kkt["stationarity_inf_norm"])),
    ("PD_alpha", float(pd["pd_alpha"])), ("PD_alpha_error", float(pd["alpha_error"])),
    ("PD_F_peak_T05_Hz", float(pd["TDS_Fpeak_T05_Hz"])),
    ("PD_R_peak_T05_Hz_s", float(pd["TDS_Rpeak_T05_Hz_s"])),
    ("support_completeness", "OPEN"), ("global_certified", "NO"),
]
pandas.DataFrame(final_metrics, columns=["metric", "value"]).to_csv(CORE / "TABLE_Q2B_final_metrics.csv", index=False)

plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "figure.dpi": 120})

# FIG_Q2B_01: the continuation stopped before 100 MW; show the direct frozen
# point separately to avoid implying a completed optimum path.
fig, ax = plt.subplots(figsize=(8.2, 4.8))
cc = continuation[continuation.disturbance_MW < 100]
ax.plot(cc.disturbance_MW, cc.cost, "o-", lw=2, color="#245b8f",
        label="CORE branch hold / continuation checkpoints")
ax.scatter([100], [J], s=95, marker="D", color="#d45a36",
           label="Direct 100 MW candidate (KKT open)", zorder=4)
ax.annotate(f"{J:.2f} MW\nnot locally certified", (100, J), xytext=(83, 445),
            textcoords="data", arrowprops={"arrowstyle": "->"}, ha="center")
ax.axvline(64.84375, ls="--", lw=1, color="#888888")
ax.text(68, 409, "continuation blocked\nnear 65.7 MW", color="#555555", fontsize=9)
ax.set(xlabel="Sustained load step at bus 16 (MW)", ylabel="Recorded retained SG (MW)",
       title="Q2B: recorded branch and frozen 100 MW candidate")
ax.set_ylim(394, 448)
fig.text(.5,.015,"Direct point is a separate all-SG-adjacent correction; no J*(d) curve or global claim.",
         ha="center", fontsize=8, color="#555555")
ax.grid(alpha=.25); ax.legend(loc="upper left", frameon=True); fig.tight_layout(rect=(0,.04,1,1))
fig.savefig(FIG / "FIG_Q2B_01_core_vs_secure.png", dpi=200); plt.close(fig)

# FIG_Q2B_02: estimator-window sensitivity at the frozen point.
fig, axs = plt.subplots(1, 2, figsize=(9.0, 4.3))
T = windows_an.window_s.to_numpy()
for ax, field, label in ((axs[0], "F_peak_Hz", "Peak |Δf_T| (Hz)"),
                         (axs[1], "R_peak_Hz_s", "Peak |R_T| (Hz/s)")):
    ax.plot(T, windows_an[field], "o-", label="Analytical", lw=2)
    ax.plot(T, windows_pd[field], "s--", label="PowerDynamics", lw=1.8)
    ax.axhline(.5, color="#b04436", ls=":", label="Project limit 0.5" if ax is axs[0] else None)
    ax.set(xlabel="Causal phase window T (s)", ylabel=label)
    ax.grid(alpha=.25)
axs[0].legend(frameon=True, fontsize=8)
fig.suptitle("Finite-window sensitivity at frozen best-found candidate")
fig.tight_layout(); fig.savefig(FIG / "FIG_Q2B_02_window_sensitivity.png", dpi=200); plt.close(fig)

# FIG_Q2B_03: fractional dispatch per generator bus.
fig, ax = plt.subplots(figsize=(8.4, 4.4))
x = np.arange(len(candidate_bus)); w = .72
ax.bar(x, candidate_bus.retained_SG_MW, w, color="#355f83", label="Retained SG MW")
ax.set_xticks(x, candidate_bus.bus.astype(str)); ax.set_xlabel("Generator bus")
ax.set_ylabel("Retained SG dispatch (MW)")
ax2 = ax.twinx(); ax2.plot(x, candidate_bus.rho, "o-", color="#d46a3b", label="ρ (GFL fraction)")
ax2.set_ylim(0, 1.02); ax2.set_ylabel("Converted fraction ρ")
ax.set_title(f"Frozen candidate: {J:.3f} MW SG retained; {GFL:.3f} MW GFL ({fraction:.3%})")
ax.grid(axis="y", alpha=.25)
lines, labels = ax.get_legend_handles_labels(); lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines+lines2, labels+labels2, loc="upper left", frameon=True)
fig.tight_layout(); fig.savefig(FIG / "FIG_Q2B_03_candidate_dispatch.png", dpi=200); plt.close(fig)

# FIG_Q2B_04: normalized security margins. Beta is highlighted because its
# observed positive margin is not interval-certified.
cm = pandas.DataFrame(constraint_rows, columns=["constraint", "observed_value", "limit", "normalized_slack", "status"])
fig, ax = plt.subplots(figsize=(8.6, 4.7))
colors = ["#477a5b" if "INTERVAL_INCOMPLETE" not in s else "#c77724" for s in cm.status]
y = np.arange(len(cm))
ax.barh(y, cm.normalized_slack, color=colors)
ax.set_yticks(y, cm.constraint); ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlabel("Normalized observed slack (log scale)")
ax.set_title("Frozen candidate constraints: observed margins")
for yy, value, status in zip(y, cm.normalized_slack, cm.status):
    ax.text(value*1.08, yy, f"{value:.2e}" + ("  (interval open)" if "INCOMPLETE" in status else ""), va="center", fontsize=8)
ax.grid(axis="x", which="both", alpha=.25); fig.tight_layout()
fig.savefig(FIG / "FIG_Q2B_04_constraint_margins.png", dpi=200); plt.close(fig)

# FIG_Q2B_05: same causal T=0.5 s operator, analytic and nonlinear PD.
pd_t = pd_trace.event_time_s.to_numpy()
phcols = [f"phase_bus{b}_rad" for b in range(30, 40)]
ph = pd_trace[phcols].to_numpy().T
dt = float(np.median(np.diff(pd_t))); lag = int(round(.5/dt)); zero = int(np.argmin(abs(pd_t)))
f_pd = np.zeros_like(ph)
for j in range(zero, ph.shape[1]):
    jm = j-lag
    theta_m = ph[:, jm] if jm >= zero else 0.0
    f_pd[:, j] = (ph[:, j]-theta_m)/(2*math.pi*.5)
analytic_t = analytic_trace.event_time_s.to_numpy()
analytic39 = analytic_trace["fGrid_bus39_T05_Hz"].to_numpy()
pd_bus39 = f_pd[9]
pd_max = np.max(np.abs(f_pd), axis=0)
an_max = analytic_trace.max_abs_fGrid_T05_Hz.to_numpy()
fig, axs = plt.subplots(2, 1, figsize=(8.5, 6.3), sharex=True)
axs[0].plot(analytic_t, analytic39, label="Analytical bus 39", lw=1.8)
axs[0].plot(pd_t, pd_bus39, "--", label="PowerDynamics bus 39", lw=1.5)
axs[0].set_ylabel("Δf_T (Hz)"); axs[0].legend(fontsize=8); axs[0].grid(alpha=.25)
axs[1].plot(analytic_t, an_max, label="Analytical max over buses 30–39", lw=1.8)
axs[1].plot(pd_t, pd_max, "--", label="PowerDynamics max over buses 30–39", lw=1.5)
axs[1].axhline(.5, color="#b04436", ls=":", label="0.5 Hz project limit")
axs[1].set(xlabel="Time from load step at bus 16 (s)", ylabel="max |Δf_T| (Hz)")
axs[1].legend(fontsize=8); axs[1].grid(alpha=.25); axs[1].set_xlim(-.5, 30)
fig.suptitle("Frozen 100 MW step, causal T=0.5 s bus-phase measurement")
fig.tight_layout(); fig.savefig(FIG / "FIG_Q2B_05_analytic_vs_PD.png", dpi=200); plt.close(fig)

# FIG_Q2B_06: local sensitivities per MW, normalized row-wise only for display.
fields = [("a_alpha_per_MW", "Modal"), ("a_beta0_per_MW", "β at ω=0"),
          ("a_Finf_per_MW", "Steady frequency"), ("a_Fpeak_per_MW", "Peak frequency"),
          ("a_R05_per_MW", "RoCoF T=0.5")]
mat = np.vstack([sensitivity[f].to_numpy() for f, _ in fields])
scale = np.max(np.abs(mat), axis=1); scale[scale == 0] = 1
normmat = mat / scale[:, None]
fig, ax = plt.subplots(figsize=(9.5, 4.4))
im = ax.imshow(normmat, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
ax.set_xticks(np.arange(10), sensitivity.bus.astype(str)); ax.set_yticks(np.arange(5), [v for _,v in fields])
for i in range(5):
    for j in range(10): ax.text(j, i, f"{mat[i,j]:.1e}", ha="center", va="center", fontsize=6.8)
ax.set_xlabel("Generator bus; local derivatives per retained fraction and MW")
ax.set_title("Local marginal sensitivities only; no KKT authority certificate")
fig.colorbar(im, ax=ax, label="Row-normalized direction (display only)")
fig.tight_layout(); fig.savefig(FIG / "FIG_Q2B_06_candidate_sensitivities.png", dpi=200); plt.close(fig)

# FIG_Q2B_07: no self-energy percentages are fabricated when alpha is slack.
fig, ax = plt.subplots(figsize=(8.2, 4.5))
ax.bar(sensitivity.bus.astype(str), sensitivity.d_alpha_d_epsilon_s_inv, color="#6987a3")
ax.axhline(0, color="#333333", lw=.8)
ax.set(xlabel="Generator bus", ylabel="dα/dε (s⁻¹)",
       title="BND decomposition gate: modal constraint is not active")
ax.text(.02,.96,
        f"α={alpha:.8f} s⁻¹; required α≤−0.05 s⁻¹\nNo modal multiplier; direct/self-energy split not reported.",
        transform=ax.transAxes, va="top", bbox={"boxstyle":"round", "fc":"white", "ec":"#aaaaaa"})
ax.grid(axis="y", alpha=.25); fig.tight_layout()
fig.savefig(FIG / "FIG_Q2B_07_BND_gate.png", dpi=200); plt.close(fig)

sha = hashlib.sha256((CORE / "Z_Q2B_SECURE_T05_FINAL.toml").read_bytes()).hexdigest()
rho_s = ", ".join(f"{v:.6f}" for v in cand["rho"])
kp_s = ", ".join(f"{v:.3f}" for v in cand["Kp"])
ki_s = ", ".join(f"{v:.3f}" for v in cand["Ki"])
window_line = "; ".join(f"T={r.window_s:g}s: R={r.R_peak_Hz_s:.6f}" for r in windows_an.itertuples())
oos_lines = "\n".join(
    "| bus {}, {} MW | {:.6f} | {:.6f} | {:.6f} | {} |".format(
        int(r["bus"]), float(r["delta_MW"]), float(r["F_peak_Hz"]),
        float(r["F_inf_Hz"]), float(r["R_peak_Hz_s"]), bool(r["pass"]))
    for r in oos.to_dict(orient="records"))
report = f"""# Experiment Q2B — result report

## Status

**EXP_Q2B_STATUS: FAIL_OPTIMIZATION** (best feasible point found, not a locally certified optimum).

The PD-exact frozen model and component P/Q contract reproduce. The frozen 100 MW candidate passes its tested PowerDynamics spectral and T=0.5 s event limits. The design search did not produce a valid KKT certificate or a completed interval certificate for the robustness radius, so no secure optimum or global optimum is claimed.

## Frozen inputs and reproducibility

- Branch: `research/expQ2B-secure-optimum`; no push or commit.
- Model SHA: `{cand['model_sha']}`.
- Julia 1.11.9, PowerDynamics {pd['versions']['PowerDynamics']}, NetworkDynamics {pd['versions']['NetworkDynamics']}.
- Candidate SHA-256: `{sha}` ([frozen TOML](CORE/Z_Q2B_SECURE_T05_FINAL.toml)). No tuning occurred after the freeze.
- The ExpQ2 regression table passes 9/9 checks; it also verifies the inherited 14/14 baseline gate, trim contract, ExpN α, continuity/PLL DC findings, and NONAFFINE steady-frequency classification.
- Loads at buses 31 and 39 were kept separate from the generator component dispatch.

## Search and candidate

The custom normalized active-set SQP used the ten-bus mixed architecture, seeded adjacent to all-SG. The d=0 checkpoint was 402.063544 MW, but was already marked `BEST_FOUND_KKT_OPEN`. At 100 MW, that d=0 point has Fpeak(T=0.5)={float(core_d0_windows.loc[core_d0_windows.window_s==.5,'F_peak_Hz'].iloc[0]):.6f} Hz, so it is event-infeasible. The held branch was feasible through 64.84375 MW; the predictor was blocked near 65.7 MW. A direct 100 MW correction from the all-SG-adjacent seed reached a feasible point. Two continuation/restoration passes reduced its cost to **{J:.6f} MW retained SG** and **{GFL:.6f} MW GFL** ({fraction:.3%} of initialized generator dispatch). The retained support is buses 30–39.

The direct candidate's normalized constraints are `{[round(float(x), 8) for x in kkt['constraint_values_scaled']]}` in order modal, robust β, steady frequency, peak frequency. Only the sampled robust row is within the 2e-4 active tolerance. The least-squares stationarity residual is **{kkt['stationarity_inf_norm']:.6f}**; LICQ rank is {kkt['LICQ_rank']} of 1; SOSC is not tested. The solver ended `NO_FEASIBLE_DESCENT`, so this is not a local optimum certificate. Only the all-SG-adjacent seed was optimized; support enumeration and multistart branch completeness remain open.

| Bus | ε retained SG | ρ GFL | Kp | Ki | Retained SG MW |
|---:|---:|---:|---:|---:|---:|
""" + "\n".join(f"| {int(r.bus)} | {r.epsilon:.6f} | {r.rho:.6f} | {r.Kp:.4f} | {r.Ki:.3f} | {r.retained_SG_MW:.4f} |" for r in candidate_bus.itertuples()) + f"""

### Frozen candidate metrics

- α analytic: {alpha:.12f} s⁻¹ (margin to −0.05 is {(-0.05-alpha):.3e} s⁻¹).
- Full frequency-refined β-star observed: {beta_obs:.15e}; requirement: {beta_req:.15e}; observed difference: {beta_obs-beta_req:.3e}.
- Peak resolvent is at ω={robust_audit['omega_peak']:.6g} rad/s. The Float64 Lipschitz audit exhausted 50,000 nodes (`INCOMPLETE_INTERVAL_BOUND`). Its partial lower estimate is {beta_interval['beta_lower_bound_float64']:.15e}, only {beta_interval['beta_lower_bound_float64']-beta_req:.3e} above the requirement, but is **not a completed certificate**. No outward-rounded proof was performed.
- Analytical bus-grid metrics: F∞={finf:.6f} Hz; Fpeak(T=0.5)={fpeak:.6f} Hz.
- Analytic finite-window RoCoF: {window_line} Hz/s. This candidate point passes all four windows; no separate window-conditioned optimum was solved because the design itself remained uncertified.
- Robustness claim: `OBSERVED_MARGIN_ONLY; CERTIFICATE_INCOMPLETE`.

## Independent PowerDynamics validation after freeze

- Trim residual max: {pd['trim_residual_max']:.3e}; component P error {pd['max_P_error_pu']:.3e} pu; Q error {pd['max_Q_error_pu']:.3e} pu; component bounds and separate ZIP loads pass.
- 203 physical finite poles; gauge pole {pd['gauge_eigenvalue'][0]:.3e}; analytic α {pd['analytic_alpha']:.12f}, PD α {pd['pd_alpha']:.12f}; max pole mismatch {pd['maximum_relevant_pole_error']:.3e} s⁻¹.
- Nonlinear +100 MW bus16 event, causal T=0.5 s: Fpeak={pd['TDS_Fpeak_T05_Hz']:.6f} Hz, Fsteady={pd['TDS_Fsteady_T05_Hz']:.6f} Hz, Rpeak={pd['TDS_Rpeak_T05_Hz_s']:.6f} Hz/s; 2% settling to the measured final plateau at {pd['settling_to_steady_after_event_s']:.2f} s. Project limits pass in this simulation.
- The analytic large-step prediction is conservative relative to PD: Fpeak differs by {pd['TDS_Fpeak_T05_Hz']-pd['analytic_Fpeak_T05_Hz']:.6f} Hz and steady frequency by {pd['TDS_Fsteady_T05_Hz']-pd['analytic_Finf_Hz']:.6f} Hz. The 100 MW response is nonlinear; the spectral identity remains excellent.

Out-of-sample PD events (same measurement, not used as constraints):

| Event | Fpeak (Hz) | Fsteady (Hz) | Rpeak (Hz/s) | Pass |
|---|---:|---:|---:|:---:|
{oos_lines}

All four out-of-sample simulations completed with successful solver retcodes. No current-limit safety claim is made.

## Why this is not an optimum claim

The robust constraint is almost exactly binding, but its interval search did not finish and the custom SQP stopped before stationarity. The KKT stationarity residual is materially nonzero; SOSC and competing supports were not certified. The d=0 continuation branch did not reach 100 MW. Therefore the result is a validated **best-found feasible candidate** for the declared event, not `Z*_secure`, not a local KKT optimum, and not global. The lower-bound/upper-bound gap remains open; the only unconditional lower bound is 0 MW retained SG.

The modal margin is slack at the candidate and there is no modal multiplier. A BND direct/self-energy decomposition was therefore withheld; the plotted bus derivatives are sensitivity diagnostics only.

## Artifacts and runtime evidence

- Design implementation: `src/bnd_expQ2B/`.
- Scripts: `experiments/bnd_expQ2B/`.
- Tables, sensor traces, candidate freeze, and audit TOMLs: `reports/experiment_Q2B/CORE/` and `REGRESSION/`.
- Figures: `reports/experiment_Q2B/FIGURES/`.
- The adaptive robustness interval audit used 50,001 nodes and 733.416 s before returning incomplete. The optimization histories record 24 initial direct iterations, 27 extension iterations (19 accepted), and 6 restoration iterations (1 accepted); no optimizer timing was instrumented. Post-freeze validation used one full-spectrum PD rebuild plus five nonlinear event simulations.

## Figures

![Continuation and frozen candidate](FIGURES/FIG_Q2B_01_core_vs_secure.png)

![Window sensitivity](FIGURES/FIG_Q2B_02_window_sensitivity.png)

![Analytic and PD response](FIGURES/FIG_Q2B_05_analytic_vs_PD.png)
"""
(REPORT / "REPORT_EXP_Q2B.md").write_text(report, encoding="utf-8")

summary = f"""# Experiment Q2B — final summary

**Status: FAIL_OPTIMIZATION.** Best-found, PowerDynamics-validated candidate; KKT and robust interval certificates remain open.

- Retained SG: **{J:.3f} MW** across buses 30–39; GFL dispatch **{GFL:.3f} MW** ({fraction:.3%}).
- Frozen candidate SHA-256: `{sha}`.
- Analytic: α={alpha:.9f} s⁻¹; β observed={beta_obs:.12e} against {beta_req:.12e}; Fpeak(T=.5)={fpeak:.6f} Hz; F∞={finf:.6f} Hz; analytic R(T=.5)={windows_an.loc[windows_an.window_s==.5,'R_peak_Hz_s'].iloc[0]:.6f} Hz/s.
- Robust interval: `INCOMPLETE_INTERVAL_BOUND` at 50,001 nodes. Candidate has only an observed beta margin, not a completed robustness certificate.
- KKT: primal {kkt['primal_violation']:.2e}; stationarity {kkt['stationarity_inf_norm']:.3f}; complementarity about {abs(kkt['multipliers_unconstrained_LS'][0]*kkt['constraint_values_scaled'][1]):.2e}; LICQ rank {kkt['LICQ_rank']}; SOSC not tested.
- Independent PD: α mismatch {pd['alpha_error']:.2e} s⁻¹; max pole mismatch {pd['maximum_relevant_pole_error']:.2e} s⁻¹; trim residual {pd['trim_residual_max']:.2e}; nonlinear bus16/100 MW Fpeak={pd['TDS_Fpeak_T05_Hz']:.4f} Hz, Rpeak={pd['TDS_Rpeak_T05_Hz_s']:.4f} Hz/s; both project event limits pass.
- All four tested causal RoCoF windows pass at this point, and the four out-of-sample PD events pass. This does not make it a locally optimal design.
- Support/global completeness: open. No global claim. No push or commit.

See [REPORT_EXP_Q2B.md](REPORT_EXP_Q2B.md) and the bundled CSV/figure evidence.
"""
(REPORT / "FINAL_SUMMARY_EXP_Q2B.md").write_text(summary, encoding="utf-8")

print(f"FIGURES={len(list(FIG.glob('*.png')))}")
print(f"Ptotal_MW={Ptotal:.9f} J_MW={J:.9f} GFL_MW={GFL:.9f} fraction={fraction:.9%}")
print(f"candidate_sha={sha}")
