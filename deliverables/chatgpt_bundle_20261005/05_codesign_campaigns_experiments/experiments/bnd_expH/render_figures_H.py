"""Render the ExpH report figures from its frozen CSV evidence."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "reports" / "experiment_H" / "tables"
FIGS = ROOT / "reports" / "experiment_H" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 200, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.22,
})


def table(name):
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGS / name, bbox_inches="tight")
    plt.close(fig)


# H01: the quotient leaves one numerically near-zero physical root for every
# deterministic gain pattern; display the residual scale rather than claiming
# a nonzero physical eigenvalue.
d = table("TABLE_H01_gauge_quotient.csv")
z = d["physical_nearzero"].astype(str).str.replace("im", "j", regex=False)
z = z.str.replace(" ", "", regex=False).str.replace("+0.0j", "", regex=False)
z = z.str.replace("+0im", "", regex=False)
vals = np.array([abs(complex(v)) for v in z])
fig, ax = plt.subplots(figsize=(7.0, 3.3))
ax.semilogy(d["pattern"], np.maximum(vals, 1e-16), "o", color="#176b87")
ax.set_ylabel(r"$|\lambda_{phys}(0,K)|$ (s$^{-1}$)")
ax.set_title("H01 · Quotient physical root remains at numerical zero")
ax.tick_params(axis="x", rotation=25)
save(fig, "FIG_H01_quotient_near_zero_spectrum.png")

# H02: small-singular-value scaling of the quotient pencil.
d = table("TABLE_H02_sigma_samples.csv").sort_values("s")
fig, ax = plt.subplots(figsize=(5.7, 3.8))
ax.loglog(d["s"], d["sigma_min"], "o-", lw=1.8, ms=4, color="#176b87")
ax.set_xlabel(r"$|s|$ (s$^{-1}$)")
ax.set_ylabel(r"$\sigma_{min}(T_q(s))$")
ax.set_title("H02 · Quotient pencil singular-value scaling")
save(fig, "FIG_H02_sigma_min_scaling.png")

# H03: exact near-zero root scaling under retained SG.
d = table("TABLE_H03_scaling_exponent.csv")
fig, ax = plt.subplots(figsize=(6.3, 4.1))
for bus, g in d.groupby("bus"):
    ax.loglog(g["epsilon"], g["abs_lambda"], color="#8ba5b1", alpha=0.52, lw=0.9)
g = d[d["bus"] == 36].sort_values("epsilon")
if len(g):
    ax.loglog(g["epsilon"], g["abs_lambda"], "o-", color="#d45d36", lw=2, ms=3,
              label=f"bus 36, fitted p={g['full_range_exponent'].iloc[0]:.3f}")
ax.loglog([1e-8, 1e-2], [1e-8, 1e-2], "--", color="#263238", label=r"reference $\epsilon^1$")
ax.set_xlabel(r"retained fraction $\epsilon$")
ax.set_ylabel(r"$|\lambda_{phys}|$ (s$^{-1}$)")
ax.set_title("H03 · Linear root law after gauge quotient")
ax.legend(frameon=False)
save(fig, "FIG_H03_puiseux_scaling.png")

# H04: singular spectrum of the normalized mixed gain Jacobian.
d = table("TABLE_H05_mixed_gain_jacobian.csv")
support = sorted(d["support_bus"].unique())
gains = sorted(d["gain_bus"].unique())
J = np.zeros((len(support), 2 * len(gains)))
for r in d.itertuples():
    i = support.index(r.support_bus)
    j = gains.index(r.gain_bus) + (len(gains) if r.parameter == "Ki" else 0)
    J[i, j] = r.derivative_normalized
s = np.linalg.svd(J, compute_uv=False)
fig, ax = plt.subplots(figsize=(5.7, 3.7))
ax.semilogy(np.arange(1, len(s) + 1), np.maximum(s, 1e-16), "o-", color="#176b87")
ax.set_xlabel("singular direction")
ax.set_ylabel("singular value of normalized Jacobian")
ax.set_title("H04 · Mixed authority is effectively rank one")
save(fig, "FIG_H04_mixed_authority_singular_values.png")

# H05: direct/self-energy cancellation in the exact authority decomposition.
d = table("TABLE_H07_BND_authority_decomposition.csv").sort_values("bus")
fig, ax = plt.subplots(figsize=(7.2, 4.0))
x = np.arange(len(d)); w = 0.38
ax.bar(x - w/2, d["A_direct"], width=w, label="direct", color="#d45d36")
ax.bar(x + w/2, d["A_self_energy"], width=w, label="collective self-energy", color="#176b87")
ax.scatter(x, d["A_total"], s=20, color="#222222", label="net authority", zorder=3)
ax.set_xticks(x, d["bus"].astype(str))
ax.set_xlabel("retained-SG support bus")
ax.set_ylabel("leading authority coefficient")
ax.set_yscale("symlog", linthresh=1)
ax.set_title("H05 · Large direct and self-energy terms largely cancel")
ax.legend(frameon=False, ncol=3, fontsize=8)
save(fig, "FIG_H05_direct_vs_self_energy_authority.png")

# H06: noncommutative graph basis versus active gain direction.
d = table("TABLE_H09_graph_active_subspace_alignment.csv")
angle = np.rad2deg(d["principal_angle_rad"].to_numpy())
captured = 100 * d["captured_active_subspace_energy"].to_numpy()
fig, ax = plt.subplots(figsize=(5.4, 3.5))
ax.bar(np.arange(len(angle)), angle, color="#d5a33b")
ax.set_ylabel("principal angle (degrees)")
ax.set_xticks(np.arange(len(angle)), [f"direction {i+1}" for i in range(len(angle))])
ax.set_title("H06 · Noncommutative graph basis captures active gain direction")
for i, c in enumerate(captured):
    ax.text(i, angle[i] + max(angle.max() * 0.04, 0.05), f"{c:.2f}% captured", ha="center", fontsize=8)
save(fig, "FIG_H06_graph_vs_active_gain_principal_angles.png")

# H07: predictor versus exact correction for each single-anchor branch.
d = table("TABLE_H13_exact_corrected_candidates.csv")
pred = d.groupby("bus")["predictor_retained_MW"].first()
exact = d[d["strict_margin"]].groupby("bus")["retained_sg_MW"].min()
common = pred.index.intersection(exact.index)
fig, ax = plt.subplots(figsize=(5.3, 4.2))
ax.scatter(pred.loc[common], exact.loc[common], color="#176b87", s=42)
for b in common:
    ax.annotate(str(b), (pred[b], exact[b]), xytext=(4, 3), textcoords="offset points", fontsize=8)
lim = max(float(pred.loc[common].max()), float(exact.loc[common].max())) * 1.04
ax.plot([0, lim], [0, lim], "--", color="#555555", lw=1)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("local/Jordan predictor retained MW")
ax.set_ylabel("exact feasible nominal retained MW")
ax.set_title("H07 · Exact all-pole correction changes the ranking")
save(fig, "FIG_H07_predictor_vs_exact_retention.png")

# H08-H09: direct resolvent robust frontier.
d = table("TABLE_H14_robustness_frontier.csv")
positive = d[d["beta_required"] > 0].sort_values("beta_required")
fig, ax = plt.subplots(figsize=(5.8, 3.7))
ax.semilogx(positive["beta_required"], positive["converted_GFL_MW"], "o-", color="#176b87")
ax.axhline(float(d.loc[d["beta_required"] == 0, "converted_GFL_MW"].iloc[0]), color="#888888", ls="--", label="nominal")
ax.set_xlabel(r"uncertainty bound $\beta$ (dimensionless)")
ax.set_ylabel("converted GFL capacity (MW)")
ax.set_title("H08 · Direct small-gain robustness frontier")
ax.legend(frameon=False)
save(fig, "FIG_H08_GFL_MW_vs_beta.png")

fig, ax = plt.subplots(figsize=(5.8, 3.7))
ax.plot(d["converted_GFL_MW"], d["small_gain_margin"], "o-", color="#d45d36")
ax.axhline(0, color="#333333", lw=1)
ax.set_xlabel("converted GFL capacity (MW)")
ax.set_ylabel(r"$1-\beta\sup_\omega\bar\sigma(M_\sigma)$")
ax.set_title("H09 · Certified small-gain margin")
save(fig, "FIG_H09_robust_margin_vs_GFL_MW.png")

# H10-H11: use independent nonlinear PD trajectories when available.
trajectory = TABLES / "TABLE_H21_tds_trajectories.csv"
if trajectory.exists():
    d = pd.read_csv(trajectory)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for pulse, g in d.groupby("pulse_fraction"):
        ax.plot(g["time_s"], g["retained_sg_coi_frequency_deviation_Hz"], label=f"pulse {pulse:g}")
    ax.set_xlabel("time (s)"); ax.set_ylabel("retained-SG COI frequency deviation (Hz)")
    ax.set_title("H10 · Nonlinear PowerDynamics frequency response")
    ax.legend(frameon=False)
    save(fig, "FIG_H10_frequency_transient.png")
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for pulse, g in d.groupby("pulse_fraction"):
        t = g["time_s"].to_numpy(); f = g["retained_sg_coi_frequency_deviation_Hz"].to_numpy()
        dt = np.diff(t); df = np.diff(f)
        keep = (dt > 0) & np.isfinite(dt) & np.isfinite(df)
        ax.plot((t[1:] + t[:-1])[keep] / 2, df[keep] / dt[keep], label=f"pulse {pulse:g}")
    ax.set_xlabel("time (s)"); ax.set_ylabel("COI RoCoF (Hz/s)")
    ax.set_title("H11 · Nonlinear PowerDynamics RoCoF")
    ax.legend(frameon=False)
    save(fig, "FIG_H11_rocof_transient.png")
else:
    d = table("TABLE_H17_transient_metrics.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.bar(["frequency", "RoCoF"], [d["unit_frequency_peak"], d["unit_rocof_peak"]], color=["#176b87", "#d45d36"])
    ax.set_ylabel("unit-disturbance peak")
    ax.set_title("H10 · Analytical linear transient peaks; TDS pending")
    save(fig, "FIG_H10_frequency_transient.png")
    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.bar(["RoCoF"], [d["unit_rocof_peak"]], color="#d45d36")
    ax.set_ylabel("unit-disturbance peak (Hz/s per MW)")
    ax.set_title("H11 · Analytical linear RoCoF peak; TDS pending")
    save(fig, "FIG_H11_rocof_transient.png")

print(f"Rendered 11 ExpH figures in {FIGS}")
