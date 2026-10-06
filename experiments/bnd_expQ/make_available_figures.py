from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
Q = ROOT / "reports" / "experiment_Q"
OUT = Q / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Fig Q1: compare the frozen replacement candidate with the independently
# validated all-SG feasibility reference using the frozen primary metric.
q1 = pd.read_csv(Q / "Q1A" / "TABLE_Q1A_frequency_metric_cases.csv")
q0 = pd.read_csv(Q / "Q0" / "TABLE_Q00_baseline_reproduction.csv")
q2 = pd.read_json(Q / "Q2_REFERENCE" / "Q2_ALL_SG_REFERENCE.json", typ="series")
expn = q1.loc[q1.scenario == "step_100MW"].iloc[0]
ref = pd.Series({
    "alpha": float(q2["alpha"]),
    "frequency": float(q2["PD_F_peak"]),
    "rocof": float(q2["PD_R_peak"]),
})
expn_alpha = float(q0.loc[q0.check == "analytic_alpha", "value"].iloc[0])
labels = ["ExpN nominal\n(frozen) ", "all-SG reference\n(feasible, not optimal)"]
alpha = [expn_alpha, ref["alpha"]]
freq = [float(expn.primary_peak_Hz), ref["frequency"]]
rocof = [float(expn.primary_RoCoF_Hz_s), ref["rocof"]]
fig, axes = plt.subplots(1, 3, figsize=(12, 4.1))
colors = ["#c65252", "#3a8f73"]
axes[0].bar(labels, alpha, color=colors)
axes[0].axhline(-0.05, color="#222", ls="--", lw=1, label="required: −0.05 s⁻¹")
axes[0].set_ylabel("Spectral abscissa (s⁻¹)")
axes[0].legend(frameon=False, fontsize=8)
axes[1].bar(labels, freq, color=colors)
axes[1].axhline(0.5, color="#222", ls="--", lw=1, label="project limit: 0.5 Hz")
axes[1].set_ylabel("Peak |local SG/PLL Δf| (Hz)")
axes[1].legend(frameon=False, fontsize=8)
axes[2].bar(labels, rocof, color=colors)
axes[2].axhline(0.5, color="#222", ls="--", lw=1, label="project limit: 0.5 Hz/s")
axes[2].set_ylabel("Peak |local SG/PLL RoCoF| (Hz/s)")
axes[2].legend(frameon=False, fontsize=8)
fig.suptitle("Stable poles do not guarantee the declared 100 MW frequency security")
for ax in axes:
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="y", alpha=.2)
fig.tight_layout()
fig.savefig(OUT / "FIG_Q1_stable_not_secure.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# Fig Q4: exact analytical reduced linear model against independent nonlinear
# PowerDynamics response for the frozen all-SG feasibility reference.
trace = pd.read_csv(Q / "Q2_REFERENCE" / "TABLE_Q2_all_sg_linear_vs_PD.csv")
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.plot(trace.time_s, trace.linear_max_abs_local_SG_Hz,
        label="ExpN analytical linear model", lw=1.7, color="#3c78a8")
ax.plot(trace.time_s, trace.PD_max_abs_local_SG_Hz,
        label="independent PowerDynamics TDS", lw=1.3, color="#d17a32", alpha=.9)
ax.axhline(.5, color="#222", ls="--", lw=1, label="project limit: 0.5 Hz")
ax.axvline(1.0, color="#888", ls=":", lw=1, label="+100 MW step at bus 16")
ax.set(xlabel="Time (s)", ylabel="Max |local SG rotor Δf| (Hz)",
       title="All-SG reference: linear prediction and nonlinear PowerDynamics")
ax.grid(alpha=.2)
ax.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(OUT / "FIG_Q4_all_sg_linear_vs_PD.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"Wrote figures to {OUT}")
