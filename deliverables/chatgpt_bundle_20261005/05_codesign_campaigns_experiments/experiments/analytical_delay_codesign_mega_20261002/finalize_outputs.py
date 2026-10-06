"""Build claim-limited tables and figures from executed Mega experiment outputs."""
from pathlib import Path
import json
import shutil
import tomllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
def csv(name):
    return pd.read_csv(HERE / name)

audit = csv("M1_CANDIDATE_AUDIT.csv")
best = audit[audit.fully_validated].sort_values("GFL_MW").iloc[-1]
design = tomllib.loads((HERE / "M1_ZERO_DELAY_DESIGN.toml").read_text(encoding="utf-8"))
historical = tomllib.loads((HERE / "m1_validations/eta_100/Q0_RESULT.toml").read_text(encoding="utf-8"))
q4 = csv("TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv")
q4.insert(0, "status", "EXACT_CONDITIONAL_TARGET_ROOT_NOT_GLOBAL_BOUNDARY")
q4.to_csv(HERE / "M5_CLOSED_FORM_PI_VALIDATION.csv", index=False)

count_parts = []
for name in ("M3_TRACE_INTEGRAL_ATTEMPT.csv", "M3_TRACE_INTEGRAL_M1_ZERO_DELAY_DESIGN.csv"):
    p = HERE / name
    if p.exists():
        frame = pd.read_csv(p)
        if "design_id" not in frame:
            frame.insert(0, "design_id", "seed_uniform_875")
        count_parts.append(frame)
counts = pd.concat(count_parts, ignore_index=True).sort_values(["design_id", "tau_ms"])
counts.to_csv(HERE / "M3_DDE_ROOT_COUNTS.csv", index=False)

root_parts = []
for name, did in (("M3_REFINED_ROOTS.csv", "seed_uniform_875"),
                  ("M3_REFINED_ROOTS_M1_ZERO_DELAY_DESIGN.csv", "M1_ZERO_DELAY_DESIGN")):
    p = HERE / name
    if p.exists():
        roots = pd.read_csv(p)
        if "design_id" not in roots:
            roots.insert(0, "design_id", did)
        root_parts.append(roots)
roots = pd.concat(root_parts, ignore_index=True) if root_parts else pd.DataFrame()
if len(roots):
    roots.to_csv(HERE / "M3_ALL_REFINED_ROOTS.csv", index=False)

oracle = []
for r in counts.itertuples():
    pairs = 0 if roots.empty else len(roots[(roots.design_id == r.design_id) & (roots.tau_ms == r.tau_ms)])
    match = int(r.nearest_integer) == 2 * pairs if int(r.nearest_integer) else True
    numerical = r.status == "STABLE_NUMERICAL_INTEGER" and match
    oracle.append({
        "design_id": r.design_id, "tau_ms": r.tau_ms,
        "roots_right_of_margin": int(r.nearest_integer),
        "phase_count": r.phase_count, "trace_count": r.trace_count_real,
        "refined_positive_imag_roots": pairs,
        "root_count_refinement_match": match,
        "status": ("NUMERICALLY_VALIDATED_STABLE" if int(r.nearest_integer) == 0 else
                   "NUMERICALLY_VALIDATED_UNSTABLE") if numerical else "INDETERMINATE",
        "directed_rounding_certified": False,
        "full_nonlinear_DDE_events_validated": False,
        "note": "Exact characteristic evaluated numerically; phase and trace agree. Not an interval certificate.",
    })
pd.DataFrame(oracle).to_csv(HERE / "M3_DDE_ORACLE.csv", index=False)

blocked = {
    "M6_CONTROLLER_AUTHORITY.csv": "No full local controller-authority prediction validated against complete DDE optimum.",
    "M7_REDUCED_OPTIMUM.csv": "Two-anchor theorem proved for reduced LP, but physical coefficients and full-model lift not validated.",
    "M8_ITERATION_TRACE.csv": "M1 is a frozen interpolation search, not the proposed analytical predictor/corrector.",
    "M9_UNIFORM_DELAY_FRONTIER.csv": "No full nonlinear DDE event solver or delay-specific co-design frontier.",
    "M10_HETEROGENEOUS_FRONTIER.csv": "No full delayed replacement frontier; placement effect cannot be inferred from fixed-design poles.",
    "M11_GRAPH_DAMPING_VALIDATION.csv": "Zero-frequency Schur block singular after gauge deflation; finite-frequency port model not validated.",
    "M12_SPATIAL_DELAY_METRICS.csv": "Frozen patterns exist; no full frontier labels for predictor validation.",
    "M13_DELAY_SHADOW_PRICES.csv": "No differentiable validated delay-specific optimum or active-set multipliers.",
    "M17_OPTIMALITY_BOUNDS.csv": "No rigorous upper replacement bound; best M1 witness is only a feasible lower bound at zero delay.",
}
for file, reason in blocked.items():
    pd.DataFrame([{"status": "BLOCKED", "reason": reason}]).to_csv(HERE / file, index=False)

pd.DataFrame([
    {"design": "uniform_875_seed", "GFL_percent": 87.5,
     "GFL_MW": 4727.415953731491, "retained_SG_MW": 675.3451362473559,
     "all_five_events_pass": True, "status": "REPRODUCED_BASELINE"},
    {"design": best.candidate_id, "GFL_percent": best.GFL_percent,
     "GFL_MW": best.GFL_MW, "retained_SG_MW": best.retained_SG_MW,
     "all_five_events_pass": True, "status": "BEST_FOUND_FROZEN_LINE_SEARCH"},
    {"design": "historical_joint_90_047", "GFL_percent": historical["GFL_percent"],
     "GFL_MW": historical["GFL_MW"], "retained_SG_MW": historical["retained_SG_MW"],
     "all_five_events_pass": False, "status": "FAILED_HOLDOUT_EVENTS"},
]).to_csv(HERE / "M15_BASELINE_COMPARISON.csv", index=False)
shutil.copyfile(HERE / "M1_ZERO_DELAY_EVENTS.csv", HERE / "M16_FULL_VALIDATION.csv")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "savefig.dpi": 200})
fig, ax = plt.subplots(figsize=(8, 4.7))
colors = np.where(audit.fully_validated, "#087f5b", "#c92a2a")
ax.scatter(audit.GFL_percent, audit.min_actuator_slack, c=colors, s=75, zorder=3)
for r in audit.itertuples():
    ax.annotate(r.candidate_id.replace("eta_", ""), (r.GFL_percent, r.min_actuator_slack),
                xytext=(4, 5), textcoords="offset points", fontsize=8)
ax.axhline(0.002, color="#db8b00", linestyle="--", label="frozen actuator slack limit")
ax.set(xlabel="GFL replacement (% of initial generator dispatch)",
       ylabel="Worst-event SG actuator slack", title="Zero-delay frozen search: five-event validation")
ax.grid(alpha=.22); ax.legend(loc="upper right", fontsize=8)
fig.tight_layout(); fig.savefig(HERE / "FIG_M1_VALIDATED_LINE_SEARCH.png"); plt.close(fig)

rank = csv("M2_ACTION_SPACE_RANK.csv")
good = rank[rank.status == "COMPLETED"]
fig, ax = plt.subplots(figsize=(7.2, 4.3))
ax.bar(["Full descriptor", "Physical dynamic", "Design action"], [282, 203, 30],
       color=["#495057", "#2b6cb0", "#087f5b"])
ax.set(ylabel="Dimension", title="Exact low-rank action on the IEEE-39 descriptor")
for i,v in enumerate((282,203,30)):
    ax.text(i,v+5,str(v),ha="center",weight="bold")
ax.text(.02,.96,f"21 accepted designs × 3 delay fields\nmax reconstruction {good.max_reconstruction_relative_error.max():.1e}",
        transform=ax.transAxes,va="top",fontsize=9)
ax.set_ylim(0,330);fig.tight_layout();fig.savefig(HERE / "FIG_F6_ACTION_SPACE_DIMENSION.png");plt.close(fig)

fig, ax = plt.subplots(figsize=(7.8, 4.6))
for did, group in counts.groupby("design_id"):
    label = "87.5% seed" if did == "seed_uniform_875" else "88.455% zero-delay design"
    offset = -.22 if did == "seed_uniform_875" else .22
    ax.scatter(group.tau_ms+offset,group.nearest_integer,s=85,label=label,zorder=3)
ax.axhline(0,color="black",lw=.8)
ax.set(xlabel="Uniform PLL measurement delay (ms)",ylabel="Characteristic roots with Re(s) > −0.05 s⁻¹",
       title="Exact DDE characteristic: numerical contour counts")
ax.set_xticks(sorted(counts.tau_ms.unique()));ax.grid(alpha=.22);ax.legend()
ax.text(.02,.98,"Sampled delays only; no delayed nonlinear events or replacement frontier",transform=ax.transAxes,
        va="top",fontsize=8,color="#a61e1e")
fig.tight_layout();fig.savefig(HERE / "FIG_M3_DDE_ROOT_COUNT.png");plt.close(fig)

if len(roots):
    fig, ax = plt.subplots(figsize=(7.8, 4.7))
    for did, group in roots.groupby("design_id"):
        label = "87.5% seed" if did == "seed_uniform_875" else "88.455% zero-delay design"
        ax.scatter(group.root_imag/(2*np.pi), group.root_real, label=label, s=80)
    ax.axhline(-.05,ls="--",color="#db8b00",label="required margin")
    ax.axhline(0,color="#444",lw=.8)
    ax.set(xlabel="Pole frequency (Hz)",ylabel="Real part (s⁻¹)",
           title="40 ms: refined unstable PLL cluster")
    ax.grid(alpha=.22);ax.legend();fig.tight_layout()
    fig.savefig(HERE / "FIG_M3_REFINED_UNSTABLE_ROOTS.png");plt.close(fig)

print("FINALIZED", "best_GFL_percent", best.GFL_percent,
      "M2_max_reconstruction", good.max_reconstruction_relative_error.max(),
      "M3_rows", len(counts), "M3_refined_positive_roots", len(roots))
