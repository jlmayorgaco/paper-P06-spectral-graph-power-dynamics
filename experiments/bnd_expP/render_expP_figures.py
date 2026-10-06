"""Render figures from saved ExpP result tables; does not rerun models."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_P"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 120})

# Conditional algebraic roots only. Color reports the complete-pole check.
roots = pd.read_csv(OUT / "P2" / "TABLE_P2_single_SG_roots.csv")
fig, ax = plt.subplots(figsize=(7.1, 4.5), constrained_layout=True)
ok = roots[roots["all_physical_poles_pass"]]
bad = roots[~roots["all_physical_poles_pass"]]
ax.scatter(ok["retained_SG_MW"], ok["alpha"], label="all physical poles pass",
           color="#167c80", marker="o", s=48)
ax.scatter(bad["retained_SG_MW"], bad["alpha"], label="another physical pole limits",
           color="#c44536", marker="x", s=58)
ax.axhline(-0.05, color="#333333", linestyle="--", linewidth=1, label="required alpha = −0.05 s⁻¹")
ax.set_xscale("log")
ax.set_xlabel("retained SG dispatch at conditional root (MW, log scale)")
ax.set_ylabel("complete-spectrum spectral abscissa (s⁻¹)")
ax.set_title("P2 conditional one-SG roots; full spectrum decides feasibility")
ax.grid(True, which="both", alpha=0.24)
ax.legend(frameon=False, fontsize=8)
fig.savefig(OUT / "FIG_P01_conditional_retention_roots.png", dpi=200)
plt.close(fig)

# Pointwise reciprocal singular-value observations are upper bounds on beta_star.
# The separate Lipschitz interval result supplies a provisional lower bound at
# one fixed-gain conditional point; its Float64 arithmetic is not outward-rounded.
obs = pd.read_csv(OUT / "P4" / "TABLE_P4_resolvent_observations.csv")
interval = pd.read_csv(OUT / "P4" / "TABLE_P4_interval_certificate_conditional.csv").iloc[0]
fig, ax = plt.subplots(figsize=(7.1, 4.5), constrained_layout=True)
ax.set_yscale("log")
ax.set_ylim(1e-12, 1e-5)
ax.axhline(float(obs["beta_req"].iloc[0]), color="#333333", linestyle="--",
           linewidth=1.1, label="required normalized beta")
for idx, row in obs.iterrows():
    color = "#c44536" if "nominal" in row["candidate"].lower() else "#d08c24"
    ax.scatter(idx, row["beta_upper_from_observation"], color=color, s=58)
    ax.annotate(row["candidate"], (idx, row["beta_upper_from_observation"]),
                xytext=((7, 7) if idx == 0 else (-170, -20)),
                textcoords="offset points", fontsize=8)
conditional_x = len(obs)
ax.scatter(conditional_x, float(interval["beta_lower_cert"]), marker="v", s=72,
           color="#167c80", label="conditional interval lower bound (Float64)")
ax.scatter(conditional_x, float(interval["beta_upper_observed"]), marker="o", s=58,
           facecolors="none", edgecolors="#d08c24", label="conditional sampled upper bound")
ax.set_xticks(range(len(obs) + 1),
              [f"ω={w:.3g} rad/s" for w in obs["omega_rad_s"]] + ["conditional point"])
ax.set_xlabel("sampled frequency witness or conditional point")
ax.set_ylabel("bounds on beta_star (pointwise upper / interval lower)")
ax.set_title("P4 observed norm and conditional interval bounds")
ax.grid(True, which="both", alpha=0.24)
ax.legend(frameon=False, fontsize=8)
margin = float(interval["beta_lower_cert"] - interval["beta_req"])
ax.text(0.98, 0.82, f"conditional lower − requirement = {margin:.2e}",
        transform=ax.transAxes, ha="right", fontsize=8, color="#167c80")
ax.text(0.02, 0.02, "Interval result uses Float64 without outward rounding; provisional only.",
        transform=ax.transAxes, fontsize=8, color="#555555")
fig.savefig(OUT / "FIG_P02_robust_pointwise_bounds.png", dpi=200)
plt.close(fig)

# Declared sustained step, not the separate 0.1-s small-pulse suite.
tds = pd.read_csv(OUT / "P5" / "TABLE_P5_declared_event_trajectories.csv")
event_t = float(tds["time_s"].min()) + 1.0
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 6.2), sharex=True,
                                constrained_layout=True)
ax1.plot(tds["time_s"], tds["COI_frequency_deviation_Hz"], color="#235789", linewidth=1.1)
ax1.axvline(event_t, color="#777777", linestyle=":", linewidth=1)
ax1.axhspan(-0.5, 0.5, color="#2a9d8f", alpha=0.12, label="±0.5 Hz project limit")
ax1.set_ylabel("COI frequency deviation (Hz)")
ax1.set_title("P5 independent PD TDS — sustained 100 MW bus-16 step")
ax1.grid(True, alpha=0.24)
ax1.legend(frameon=False, fontsize=8, loc="best")
ax2.plot(tds["time_s"], tds["RoCoF_Hz_s"], color="#c44536", linewidth=0.9)
ax2.axvline(event_t, color="#777777", linestyle=":", linewidth=1)
ax2.axhspan(-0.5, 0.5, color="#2a9d8f", alpha=0.12, label="±0.5 Hz/s project limit")
ax2.set_xlabel("time (s)")
ax2.set_ylabel("sampled RoCoF (Hz/s)")
ax2.grid(True, alpha=0.24)
ax2.legend(frameon=False, fontsize=8, loc="best")
fig.savefig(OUT / "FIG_P03_declared_step_tds.png", dpi=200)
plt.close(fig)

print("Rendered FIG_P01_conditional_retention_roots.png, FIG_P02_robust_pointwise_bounds.png, FIG_P03_declared_step_tds.png")
