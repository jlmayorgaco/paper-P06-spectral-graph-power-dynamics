from pathlib import Path
import json

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

csv_path = ROOT / "outputs" / "phase_braess_nep" / "phase_braess_nep_line_decomposition.csv"
status_path = ROOT / "outputs" / "phase_braess_nep" / "phase_braess_nep_status.json"

df = pd.read_csv(csv_path)
with open(status_path, "r", encoding="utf-8") as f:
    status = json.load(f)

case = df[
    (df["case"] == "mix60_no_pss")
    & (df["line_idx"] == "Line_42")
    & (df["from_bus"] == 20)
    & (df["to_bus"] == 34)
    & (df["mode_rank"] == 2)
].iloc[0]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# Figure 1: application diagnostic report for the clearest weak link.
fig = plt.figure(figsize=(7.2, 3.1))
gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1.1, 1.2], height_ratios=[1, 1])

ax0 = fig.add_subplot(gs[:, 0])
ax0.axis("off")
report_lines = [
    "Weak-link report",
    "Case: Mix60/no-PSS",
    "Action: strengthen Line 20-34 by 1%",
    "",
    "Classical view:",
    "  graph / stiffness terms nearly silent",
    "  A = -1.09e-9, B = +7.97e-9",
    "",
    "Controller-aware view:",
    "  PLL-dominated control pole is close",
    "  self-energy term C = -6.87e-3",
    "  control fraction = 0.959",
    "",
    "Decision:",
    "  do not reinforce this line alone;",
    "  retune/control-separate PLL pole or add damping."
]
ax0.text(0.02, 0.98, "\n".join(report_lines), va="top", ha="left",
         fontsize=7.5,
         bbox=dict(boxstyle="round,pad=0.35", fc="#f8fafc", ec="#334155", lw=0.8))

ax1 = fig.add_subplot(gs[0, 1:])
terms = ["A stiffness", "B retained", "C self-energy"]
vals = [case["A_stiffness_dzeta"], case["B_retained_rotation_dzeta"], case["C_self_energy_dzeta"]]
colors = ["#94a3b8", "#60a5fa", "#ef4444"]
ax1.barh(terms, vals, color=colors)
ax1.axvline(0, color="black", lw=0.8)
ax1.set_title("Why the reinforcement is weak: C dominates")
ax1.set_xlabel(r"$\partial\zeta/\partial w$ contribution")
ax1.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))

ax2 = fig.add_subplot(gs[1, 1])
dz = case["post_network_zeta"] - case["base_network_zeta"]
ax2.bar([""], [dz * 1e5], color="#dc2626")
ax2.axhline(0, color="black", lw=0.8)
ax2.set_title("Damping margin falls")
ax2.set_ylabel("")
ax2.text(0, dz * 1e5 * 0.55, r"$\Delta\zeta=-6.85{\times}10^{-5}$",
         ha="center", va="center", color="white", fontsize=7, weight="bold")
ax2.grid(alpha=0.25)

ax3 = fig.add_subplot(gs[1, 2])
ddist = case["distance_to_acc_post"] - case["distance_to_acc_base"]
ax3.bar([""], [ddist * 1e3], color="#7c3aed")
ax3.axhline(0, color="black", lw=0.8)
ax3.set_title("Control pole gets closer")
ax3.set_ylabel("")
ax3.text(0, ddist * 1e3 * 0.55, r"$\Delta d=-4.37{\times}10^{-3}$",
         ha="center", va="center", color="white", fontsize=7, weight="bold")
ax3.grid(alpha=0.25)

save(fig, "fig_application_weaklink_report")


# Figure 2: audit-level counts and mechanism taxonomy.
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
counts = {
    "all tests": status["n_rows"],
    "global\nreversals": status["n_margin_reversals"],
    "network-mode\nreversals": status["n_network_mode_reversals"],
    "strict\nresonant": status["n_resonant"],
}
axes[0].bar(counts.keys(), counts.values(), color=["#64748b", "#f97316", "#fb923c", "#dc2626"])
axes[0].set_title("Audit funnel: from line tests to strict resonance")
axes[0].set_ylabel("records")
axes[0].grid(axis="y", alpha=0.25)

class_counts = df["classification"].value_counts()
order = ["geometric_or_stiffness", "self_energy_nonresonant", "resonant", "nonreversal_or_other"]
labels = []
sizes = []
for key in order:
    if key in class_counts:
        labels.append(key.replace("_", "\n"))
        sizes.append(class_counts[key])
axes[1].bar(labels, sizes, color=["#94a3b8", "#60a5fa", "#ef4444", "#cbd5e1"][:len(labels)])
axes[1].set_title("Mechanism labels")
axes[1].set_ylabel("records")
axes[1].grid(axis="y", alpha=0.25)

save(fig, "fig_audit_mechanism_summary")


# Figure 3: method diagram, intentionally simple for conference readability.
fig, ax = plt.subplots(figsize=(7.2, 2.25))
ax.axis("off")
boxes = [
    ("Static graph\ncentrality", "misses controller poles", "#e2e8f0"),
    ("Full ANDES\nlinearization", "separate network and controls", "#dbeafe"),
    ("Schur NEP\nbridge", "keep control resolvent", "#ede9fe"),
    ("Weak-element\ncertificate", "damping loss + pole proximity", "#fee2e2"),
    ("Robustness\naction", "avoid line-only fix;\nretune or add damping", "#dcfce7"),
]
xpos = [0.03, 0.245, 0.46, 0.675, 0.88]
for i, ((title, subtitle, color), x) in enumerate(zip(boxes, xpos)):
    ax.text(x, 0.60, title, ha="center", va="center", weight="bold",
            bbox=dict(boxstyle="round,pad=0.35", fc=color, ec="#334155", lw=0.8))
    ax.text(x, 0.20, subtitle, ha="center", va="center", fontsize=7)
    if i < len(boxes) - 1:
        ax.annotate("", xy=(xpos[i + 1] - 0.075, 0.60), xytext=(x + 0.075, 0.60),
                    arrowprops=dict(arrowstyle="->", lw=1.1, color="#334155"))
ax.text(0.02, 0.92, "Application workflow: why graph-only analysis misses the hidden margin",
        ha="left", va="center", fontsize=10, weight="bold")
save(fig, "fig_application_workflow")


report = f"""# Weak-link application report

Case: Mix60/no-PSS
Line: {case['line_idx']} ({int(case['from_bus'])}-{int(case['to_bus'])})
Perturbation: 1% line strengthening

Observed effect:
- Global margin zeta: {case['base_margin_zeta']:.9f} -> {case['post_margin_zeta']:.9f}
- Tracked network-mode zeta: {case['base_network_zeta']:.9f} -> {case['post_network_zeta']:.9f}
- Distance to matched control pole: {case['distance_to_acc_base']:.6f} -> {case['distance_to_acc_post']:.6f}

Why static analysis fails:
- Retained stiffness term A: {case['A_stiffness_dzeta']:.3e}
- Retained rotation term B: {case['B_retained_rotation_dzeta']:.3e}
- Controller self-energy term C: {case['C_self_energy_dzeta']:.3e}
- Control denominator fraction: {case['control_denominator_fraction']:.3f}

Engineering interpretation:
The line is weak for this action because strengthening it pulls a network-family
mode toward a PLL-dominated condensed-control pole. A graph-only or stiffness-only
ranking does not contain that pole distance and therefore misses the hidden margin.

Recommended action:
Do not apply line-only reinforcement as the first fix. Retune the implicated PLL/control
pole, add damping, or choose a reinforcement that increases damping without reducing
the network-control pole distance.
"""
(Path(__file__).resolve().parent / "application_report_weaklink_line20_34.md").write_text(report, encoding="utf-8")

print("generated application figures and report")
