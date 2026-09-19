"""Generate compact campaign figures from recorded Gate-0 tables."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
TABLES = CAMPAIGN / "derived" / "tables"
PNG = CAMPAIGN / "figures" / "png"
PDF = CAMPAIGN / "figures" / "pdf"
SVG = CAMPAIGN / "figures" / "svg"
for directory in (PNG, PDF, SVG):
    directory.mkdir(parents=True, exist_ok=True)


def read(name: str) -> list[dict[str, str]]:
    with (TABLES / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save(fig, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(PNG / f"{stem}.png", dpi=220)
    fig.savefig(PDF / f"{stem}.pdf")
    fig.savefig(SVG / f"{stem}.svg")
    plt.close(fig)


v4 = read("TABLE_V4_PORTFOLIOS.csv")
labels = [row["metric"].replace(" [s^-1]", "") for row in v4]
values = [float(row["observed"]) for row in v4[:2]]
fig, ax = plt.subplots(figsize=(7.0, 4.0))
bars = ax.bar(["P4 nominal H4", "P4 governed H4"], values, color=["#b42318", "#1f7a5c"])
ax.axhline(0, color="#222222", linewidth=0.8)
ax.set_ylabel("alpha [s^-1]")
ax.set_title("Frozen P4 H4 policy dependence")
ax.text(0, values[0], f" {values[0]:.6f}", va="bottom", ha="center", fontsize=9)
ax.text(1, values[1], f" {values[1]:.6f}", va="top", ha="center", fontsize=9)
save(fig, "F1_p4_policy_dependence")

v9 = read("TABLE_V9_BLOCKERS.csv")
lookup = {
    row["metric"]: float(row["observed"])
    for row in v9
    if row["metric"] != "V9 same-policy global false-safe class"
}
fig, ax = plt.subplots(figsize=(7.0, 4.0))
x = [0, 1, 2, 3]
values = [lookup["V9 same-policy global correct"], lookup["V9 target 0.3-1.5 Hz correct"],
          lookup["V9 same-policy global false-safe"], lookup["V9 target 0.3-1.5 Hz false-safe"]]
bars = ax.bar(x, values, color=["#d97706", "#1f7a5c", "#b42318", "#2563eb"])
ax.set_xticks(x, ["global\ncorrect", "targeted\ncorrect", "global\nfalse-safe", "targeted\nfalse-safe"])
ax.set_ylabel("count")
ax.set_title("Corrected same-policy V9 scope split")
for bar, value in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width() / 2, value, f" {value:.0f}", ha="center", va="bottom", fontsize=9)
save(fig, "F2_v9_transfer_refutation")

tds = read("TABLE_TDS.csv")
tds_observed = [float(row["observed"]) for row in tds]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.bar(["declared agreements", "P4 H4 unstable rows"], tds_observed, color=["#2563eb", "#b42318"])
ax.set_ylabel("count")
ax.set_title("Frozen nonlinear/TDS audit")
for idx, value in enumerate(tds_observed):
    ax.text(idx, value, f" {value:.0f}", ha="center", va="bottom")
save(fig, "F3_frozen_tds_audit")

(CAMPAIGN / "figures" / "FIGURE_MANIFEST.md").write_text(
    "# Figure Manifest\n\n"
    "F1 and F2 are generated from the compact Gate-0 tables. F3 is generated\n"
    "from the frozen TDS summary. These figures are campaign artifacts; they\n"
    "do not imply that the NOT_TESTED robustness and holdout gates passed.\n",
    encoding="utf-8",
)
print("generated F1-F3 in png/pdf/svg")
