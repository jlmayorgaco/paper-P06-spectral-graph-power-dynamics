from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "PD39_blocker_atlas"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 10,
    "figure.dpi": 160,
    "savefig.dpi": 240,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def save(fig, stem):
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def read_csv(name):
    with (ROOT / "results" / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# Same selected cases under the default census path and strict Gate A path.
cases = read_csv("PD39_BLOCKER_ATLAS_MASTER_CASES.csv")
selected = [r for r in cases if r["role"] in ("selected_blocker", "matched_control")]
labels = ["H1 blocker\n35;36", "H1 control\n32;36", "H2 blocker\n37;38", "H2 control\n36;38", "H3 blocker\n30;32;33", "H3 control\n30;33;35"]
lookup = {(r["role"], r["portfolio"]): r for r in selected}
keys = [("selected_blocker", "35;36"), ("matched_control", "32;36"), ("selected_blocker", "37;38"), ("matched_control", "36;38"), ("selected_blocker", "30;32;33"), ("matched_control", "30;33;35")]
default = [float(lookup[k]["default_path_alpha"]) for k in keys]
strict = [float(lookup[k]["strict_1e10_alpha"]) for k in keys]
x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(10.2, 4.9))
w = 0.36
ax.bar(x - w / 2, default, w, color="#c65353", label="Default census initialization")
ax.bar(x + w / 2, strict, w, color="#2c6e9e", label="Gate A path at 1e-10")
ax.axhline(0, color="#222", lw=0.9)
ax.axhline(-0.05, color="#666", lw=0.9, ls="--", label="-0.05 line")
ax.set_xticks(x, labels)
ax.set_ylabel("Rightmost eigenvalue alpha [s^-1]")
ax.set_title("PD39 physicality audit: path-dependent blocker signs")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14))
ax.set_ylim(min(strict + default) - 2, max(strict + default) + 2)
for i, (a, b) in enumerate(zip(default, strict)):
    ax.text(i - w / 2, a + (0.45 if a >= 0 else -0.45), f"{a:.2f}", ha="center", va="bottom" if a >= 0 else "top", fontsize=8)
    ax.text(i + w / 2, b + (0.45 if b >= 0 else -0.45), f"{b:.2f}", ha="center", va="bottom" if b >= 0 else "top", fontsize=8)
save(fig, "F1_path_dependence_selected_cases")


# Common 1% P/Q TDS at C12. Failed blockers are explicitly shown as MaxIters.
tds = [r for r in read_csv("PD39_BLOCKER_TDS_VALIDATION.csv") if r["condition"] == "C12"]
ordered = ["H1_blocker", "H1_pred1", "H1_pred2", "H1_control", "H2_blocker", "H2_pred1", "H2_pred2", "H2_control", "H3_blocker", "H3_pred1", "H3_pred2", "H3_pred3", "H3_control"]
tds_lookup = {r["case_id"]: r for r in tds}
vals = []
colors = []
for case_id in ordered:
    r = tds_lookup[case_id]
    val = float(r["nonlinear_rate_s_inv"]) if r["nonlinear_rate_s_inv"] not in ("", "NaN") else np.nan
    vals.append(val)
    colors.append("#c65353" if "blocker" in r["role"] else "#2c6e9e")
fig, ax = plt.subplots(figsize=(10.2, 4.8))
x = np.arange(len(ordered))
plot_vals = [v if math.isfinite(v) else 0 for v in vals]
ax.bar(x, plot_vals, color=colors)
for i, v in enumerate(vals):
    if not math.isfinite(v):
        ax.text(i, 0.006, "MaxIters", ha="center", va="bottom", rotation=90, fontsize=8, color="#9a2d2d")
ax.axhline(0, color="#222", lw=0.9)
ax.set_xticks(x, [s.replace("_", "\n") for s in ordered], fontsize=8)
ax.set_ylabel("Estimated nonlinear rate [s^-1]")
ax.set_title("Common 1% P/Q TDS at C12: blockers do not complete")
ax.text(0.01, 0.98, "Red: selected blocker; blue: predecessor/control.\nBlocker integrations returned the preregistered max-iteration failure.", transform=ax.transAxes, va="top", fontsize=9)
save(fig, "F2_tds_stop_summary")


# Gate A tolerance sweep for the three selected C12 portfolios.
audit = read_csv("PD39_BLOCKER_NUMERICAL_AUDIT.csv")
targets = ["35;36", "37;38", "30;32;33"]
fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), sharey=True)
tol_order = ["1.0e-8", "1.0e-10", "1.0e-12"]
for ax, p in zip(axes, targets):
    rr = [r for r in audit if r["portfolio"] == p and r["condition"] == "C12"]
    native = {r["tolerance"]: float(r["native_alpha"]) for r in rr}
    ax.plot(np.arange(3), [native.get(t, np.nan) for t in tol_order], marker="o", color="#2c6e9e", lw=2)
    ax.axhline(0, color="#222", lw=0.8)
    ax.set_xticks(np.arange(3), ["1e-8", "1e-10", "1e-12"])
    ax.set_title(p)
    ax.set_xlabel("Equilibrium tolerance")
    ax.grid(axis="y", alpha=0.22)
axes[0].set_ylabel("Native alpha [s^-1]")
fig.suptitle("Gate A tolerance sweep: selected C12 candidates", y=1.02)
save(fig, "F3_gate_a_tolerance_sweep")

print(f"wrote {OUT}")
