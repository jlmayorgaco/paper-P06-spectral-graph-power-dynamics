from __future__ import annotations

import json
import math
import subprocess
import sys
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
POSTER = ROOT / "poster_ieee_session"
DATA = POSTER / "data"
FIG = POSTER / "figures"

BLUE = "#00843D"
DARK = "#182A20"
CYAN = "#2FAE73"
GREEN = "#00843D"
RED = "#B23A48"
ORANGE = "#D98C1F"
GRAY = "#4F6F5D"
LIGHT = "#E4F4EB"
PAPER = "#FFFFFF"


def ensure_demo_data() -> None:
    if not (DATA / "demo_case_summary.json").exists():
        subprocess.run(
            [sys.executable, str(POSTER / "scripts" / "generate_demo_and_figures.py")],
            cwd=ROOT,
            check=True,
        )


def load_case() -> dict:
    ensure_demo_data()
    with (DATA / "demo_case_summary.json").open("r", encoding="utf-8") as f:
        return json.load(f)


def save(fig: plt.Figure, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=360, bbox_inches="tight")
    plt.close(fig)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": DARK,
            "axes.labelcolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
            "text.color": DARK,
            "axes.titleweight": "bold",
            "axes.titlesize": 16,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.facecolor": PAPER,
            "axes.facecolor": PAPER,
        }
    )


def line_name(row: dict) -> str:
    return f"{row['from_bus']}-{row['to_bus']}"


def load_rankings() -> list[dict]:
    with (DATA / "demo_line_rankings.csv").open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    int_keys = {
        "edge_index",
        "from_bus",
        "to_bus",
        "frequency_rank",
        "damping_rank",
    }
    bool_keys = {
        "is_frequency_top",
        "is_damping_top",
        "is_reversal",
        "second_order_direction_correct",
    }
    for row in rows:
        for key, value in list(row.items()):
            if key in int_keys:
                row[key] = int(value)
            elif key in bool_keys:
                row[key] = value == "True"
            else:
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    pass
    return rows


def plot_witness_bar_ranking(case: dict) -> None:
    rows = load_rankings()
    top_freq = case["selected_reversal_frequency_top"]
    top_damp = case["selected_damping_aware_top"]

    def color_for(row: dict) -> str:
        if row["edge_index"] == top_freq["edge_index"]:
            return RED
        if row["edge_index"] == top_damp["edge_index"]:
            return GREEN
        if row["delta_zeta_full_fd_10pct"] < 0:
            return "#B8C2CC"
        return "#9BB8AA"

    frequency_rows = sorted(rows, key=lambda r: r["delta_nu_fd_10pct"], reverse=True)
    damping_rows = sorted(rows, key=lambda r: r["delta_zeta_full_fd_10pct"], reverse=True)

    fig, axes = plt.subplots(1, 2, figsize=(18.5, 13.4), gridspec_kw={"wspace": 0.38})
    fig.patch.set_facecolor(PAPER)

    panels = [
        (
            axes[0],
            frequency_rows,
            "Frequency ranking",
            r"critical stiffness movement $\Delta\nu_c$",
            "delta_nu_fd_10pct",
            False,
        ),
        (
            axes[1],
            damping_rows,
            "Damping-aware ranking",
            r"verified margin movement $\Delta\zeta_{\min}$",
            "delta_zeta_full_fd_10pct",
            True,
        ),
    ]

    for ax, panel_rows, title, xlabel, value_key, show_zero in panels:
        labels = [line_name(row) for row in panel_rows]
        values = np.array([row[value_key] for row in panel_rows], dtype=float)
        y = np.arange(len(panel_rows))
        colors = [color_for(row) for row in panel_rows]

        ax.barh(y, values, color=colors, edgecolor="white", linewidth=2.2, height=0.72)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=23, fontweight="bold")
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=29, fontweight="bold", pad=18)
        ax.set_xlabel(xlabel, fontsize=24, labelpad=12)
        ax.tick_params(axis="x", labelsize=20)
        ax.grid(axis="x", alpha=0.22, linewidth=1.4)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.8)
        if show_zero:
            ax.axvline(0.0, color=DARK, linewidth=2.0)

        span = max(values.max() - values.min(), np.max(np.abs(values)), 1e-6)
        offset = 0.025 * span
        for yi, row, value in zip(y, panel_rows, values):
            if value_key == "delta_nu_fd_10pct":
                label = f"{value:+.3f}"
            else:
                label = f"{value:+.4f}"
            ha = "left" if value >= 0 else "right"
            x_text = value + offset if value >= 0 else value - offset
            ax.text(
                x_text,
                yi,
                label,
                va="center",
                ha=ha,
                fontsize=20,
                fontweight="bold" if row["edge_index"] in {top_freq["edge_index"], top_damp["edge_index"]} else "normal",
                color=DARK,
            )

    axes[0].text(
        top_freq["delta_nu_fd_10pct"] * 0.56,
        [line_name(row) for row in frequency_rows].index(line_name(top_freq)),
        "HARMFUL\nfull-QEP damping drops",
        va="center",
        ha="center",
        fontsize=20,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", fc=RED, ec=RED),
    )
    axes[1].text(
        top_damp["delta_zeta_full_fd_10pct"] * 0.52,
        [line_name(row) for row in damping_rows].index(line_name(top_damp)),
        "CHOSEN\nmargin improves",
        va="center",
        ha="center",
        fontsize=20,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", fc=GREEN, ec=GREEN),
    )

    axes[1].set_xlim(
        min(row["delta_zeta_full_fd_10pct"] for row in rows) * 1.35,
        max(row["delta_zeta_full_fd_10pct"] for row in rows) * 1.55,
    )
    axes[0].set_xlim(0, max(row["delta_nu_fd_10pct"] for row in rows) * 1.20)

    fig.suptitle(
        "Same candidate lines, different decision: stiffness ranking is not a damping-margin certificate",
        fontsize=31,
        fontweight="bold",
        color=DARK,
        y=1.02,
    )
    fig.text(
        0.50,
        -0.005,
        (
            f"Master seed {case['master_seed']} | Case seed {case['case_seed']} | "
            f"{case['reversal_count']}/{case['n_lines_tested']} tested lines raise stiffness but lower damping margin"
        ),
        ha="center",
        fontsize=19,
        color=GRAY,
    )
    save(fig, "fig_witness_bar_ranking")


def plot_reversal_story(case: dict) -> None:
    weights = np.array(case["weight_matrix"], dtype=float)
    inertia = np.array(case["inertia_M"], dtype=float)
    damping = np.array(case["damping_per_inertia_D_over_M"], dtype=float)
    rows = sorted(case["edges"], key=lambda r: r["edge_index"])

    ranking_rows = load_rankings()

    top_freq = case["selected_reversal_frequency_top"]
    top_damp = case["selected_damping_aware_top"]

    fig = plt.figure(figsize=(18.2, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.28, 1.10], wspace=0.48)

    ax_net = fig.add_subplot(gs[0, 0])
    graph = nx.Graph()
    for i in range(weights.shape[0]):
        graph.add_node(i + 1)
    for row in rows:
        graph.add_edge(row["from_bus"], row["to_bus"], weight=row["base_weight"])

    pos = nx.spring_layout(graph, seed=17, weight="weight", iterations=300)
    edge_colors = []
    edge_widths = []
    for u, v in graph.edges():
        pair = {u, v}
        if pair == {top_freq["from_bus"], top_freq["to_bus"]}:
            edge_colors.append(RED)
            edge_widths.append(5.5)
        elif pair == {top_damp["from_bus"], top_damp["to_bus"]}:
            edge_colors.append(GREEN)
            edge_widths.append(5.5)
        else:
            edge_colors.append("#B7C1CE")
            edge_widths.append(1.8)

    nx.draw_networkx_edges(graph, pos, ax=ax_net, edge_color=edge_colors, width=edge_widths)
    node_sizes = 900 + 900 * (inertia - inertia.min()) / (inertia.max() - inertia.min())
    nodes = nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax_net,
        node_size=node_sizes,
        node_color=damping,
        cmap="Blues",
        edgecolors=DARK,
        linewidths=1.8,
    )
    nx.draw_networkx_labels(graph, pos, ax=ax_net, font_color="white", font_weight="bold")
    ax_net.set_title("A. Six-bus witness network", loc="left")
    ax_net.text(
        0.02,
        -0.08,
        "Node color: local damping $D_i/M_i$   Node size: inertia $M_i$",
        transform=ax_net.transAxes,
        fontsize=9.5,
        color=GRAY,
    )
    ax_net.text(
        0.02,
        0.98,
        "red = frequency pick\n"
        "green = damping pick",
        transform=ax_net.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#D2D9E2"),
    )
    ax_net.set_axis_off()
    cbar = fig.colorbar(nodes, ax=ax_net, shrink=0.62, pad=0.02)
    cbar.set_label(r"$D_i/M_i$", fontsize=10)

    ax_scatter = fig.add_subplot(gs[0, 1])
    x = np.array([r["dnu_dw_frequency"] for r in ranking_rows])
    y = np.array([r["delta_zeta_full_fd_10pct"] for r in ranking_rows])
    colors = np.where(y < 0, RED, GREEN)
    ax_scatter.axhline(0.0, color=DARK, lw=1.0)
    ax_scatter.scatter(x, y, s=110, c=colors, edgecolors="white", linewidths=1.2, alpha=0.95)
    ax_scatter.scatter(
        [top_freq["dnu_dw_frequency"]],
        [top_freq["delta_zeta_full_fd_10pct"]],
        s=260,
        c=RED,
        marker="X",
        edgecolors=DARK,
        linewidths=1.2,
        label=f"frequency top: {line_name(top_freq)}",
    )
    ax_scatter.scatter(
        [top_damp["dnu_dw_frequency"]],
        [top_damp["delta_zeta_full_fd_10pct"]],
        s=240,
        c=GREEN,
        marker="*",
        edgecolors=DARK,
        linewidths=1.0,
        label=f"damping top: {line_name(top_damp)}",
    )
    for r in ranking_rows:
        if r["edge_index"] in {top_freq["edge_index"], top_damp["edge_index"]}:
            ax_scatter.annotate(
                line_name(r),
                (r["dnu_dw_frequency"], r["delta_zeta_full_fd_10pct"]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=11,
                fontweight="bold",
            )
    ax_scatter.set_xlabel(r"frequency-only score $\partial\nu_c/\partial w$")
    ax_scatter.set_ylabel(r"verified full-QEP movement $\Delta\zeta_{\min}$")
    ax_scatter.set_title("B. Frequency can rank harmful reinforcements first", loc="left")
    ax_scatter.grid(True, alpha=0.22)
    ax_scatter.legend(frameon=False, loc="lower right", fontsize=9.5)
    ax_scatter.text(
        0.02,
        0.95,
        f"{case['reversal_count']}/{case['n_lines_tested']} tested lines raise $\\nu_c$ but lower $\\zeta_{{\\min}}$",
        transform=ax_scatter.transAxes,
        va="top",
        fontsize=12,
        color=RED,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="#D2D9E2", alpha=0.92),
    )

    ax_cards = fig.add_subplot(gs[0, 2])
    ax_cards.set_axis_off()
    ax_cards.set_title("C. What the method catches", loc="left", pad=20)
    z0 = case["base"]["zeta_full"]
    z1 = top_freq["post_zeta_full"]
    z2 = top_damp["post_zeta_full"]
    vals = [z0, z1, z2]
    labels = ["base", "frequency pick\n1-3", "damping pick\n2-5"]
    xs = np.array([0.15, 0.5, 0.85])
    y0 = 0.50
    ax_cards.plot([xs[0], xs[1]], [y0, y0 - 0.16], color=RED, lw=4, marker="o", ms=9)
    ax_cards.plot([xs[0], xs[2]], [y0, y0 + 0.12], color=GREEN, lw=4, marker="o", ms=9)
    for xi, val, lab, dy in zip(xs, vals, labels, [0, -0.16, 0.12]):
        ax_cards.text(xi, y0 + dy + 0.055, f"{val:.5f}", ha="center", va="bottom", fontsize=17, fontweight="bold")
        ax_cards.text(xi, y0 + dy - 0.065, lab, ha="center", va="top", fontsize=11)
    ax_cards.text(
        0.5,
        0.85,
        r"Full-QEP reference confirms the reversal",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color=DARK,
    )
    ax_cards.text(
        0.5,
        0.18,
        (
            "Second order predicts the harmful direction\n"
            f"$\\Delta\\widehat{{\\zeta}}_2={top_freq['delta_zeta_second_order_fd_10pct']:+.4f}$ "
            f"vs. full-QEP $\\Delta\\zeta={top_freq['delta_zeta_full_fd_10pct']:+.4f}$"
        ),
        ha="center",
        va="center",
        fontsize=12.5,
        bbox=dict(boxstyle="round,pad=0.55", fc=LIGHT, ec="#D2D9E2"),
    )
    ax_cards.set_xlim(0, 1)
    ax_cards.set_ylim(0, 1)

    fig.suptitle(
        "Line reinforcement reversal: a cheap frequency heuristic chooses a line that worsens damping",
        fontsize=21,
        fontweight="bold",
        color=DARK,
        y=1.03,
    )
    save(fig, "fig_reversal_story")


def plot_method_pipeline() -> None:
    fig, ax = plt.subplots(figsize=(15.2, 4.1))
    ax.set_axis_off()
    steps = [
        ("1", "Normalize", r"$\widetilde L=M^{-1/2}LM^{-1/2}$" "\n" r"$\widetilde D=M^{-1/2}DM^{-1/2}$"),
        ("2", "Graph modes", r"$\widetilde Lq_k=\nu_k q_k$" "\n" r"$\Gamma=Q^T\widetilde DQ$"),
        ("3", "Diagonal screen", r"$\widehat\zeta_k=\Gamma_{kk}/(2\sqrt{\nu_k})$" "\n" "fast, not a certificate"),
        ("4", "Rational correction", r"$\Delta s_c\propto [E\,H(\Lambda;s_0)\,E]_{cc}$" "\n" "captures modal damping coupling"),
        ("5", "Fallback", "if modes cluster or correction is large\nsolve reduced QEP / full QEP"),
    ]
    xs = np.linspace(0.08, 0.92, len(steps))
    for idx, (num, title, body) in enumerate(steps):
        x = xs[idx]
        ax.add_patch(
            plt.Rectangle(
                (x - 0.080, 0.25),
                0.16,
                0.50,
                facecolor=LIGHT if idx != 3 else "#E5F3FF",
                edgecolor=BLUE if idx == 3 else "#B7C1CE",
                linewidth=2.2,
                transform=ax.transAxes,
            )
        )
        ax.add_patch(
            plt.Circle((x - 0.067, 0.69), 0.028, transform=ax.transAxes, facecolor=BLUE, edgecolor="white", lw=1.2)
        )
        ax.text(x - 0.067, 0.69, num, transform=ax.transAxes, ha="center", va="center", color="white", fontweight="bold")
        ax.text(x, 0.60, title, transform=ax.transAxes, ha="center", va="center", fontsize=13.0, fontweight="bold")
        ax.text(x, 0.43, body, transform=ax.transAxes, ha="center", va="center", fontsize=9.8)
        if idx < len(steps) - 1:
            ax.annotate(
                "",
                xy=(xs[idx + 1] - 0.098, 0.50),
                xytext=(x + 0.090, 0.50),
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", lw=2.2, color=GRAY),
            )
    ax.text(
        0.5,
        0.08,
        "Core idea: network frequency alone is diagonal in the graph basis; IBR damping enters through off-diagonal modal coupling.",
        transform=ax.transAxes,
        ha="center",
        fontsize=13.5,
        fontweight="bold",
        color=DARK,
    )
    save(fig, "fig_method_pipeline")


def plot_evidence_compact(case: dict) -> None:
    with (ROOT / "outputs" / "ieee39_rational_filter" / "surrogate_summary.json").open("r", encoding="utf-8") as f:
        surrogate = json.load(f)

    fig = plt.figure(figsize=(16.2, 5.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.95, 1.05], wspace=0.35)

    ax0 = fig.add_subplot(gs[0, 0])
    regimes = ["weak", "moderate", "strong", "clustered"]
    diag = np.array([0.038, 2.53, 5.16, 308.0])
    second = np.array([0.007, 0.65, 1.46, 74.0])
    x = np.arange(len(regimes))
    width = 0.34
    ax0.bar(x - width / 2, diag, width=width, color="#87AFCB", label="diagonal")
    ax0.bar(x + width / 2, second, width=width, color=BLUE, label="second order")
    ax0.set_yscale("log")
    ax0.set_xticks(x)
    ax0.set_xticklabels(regimes, rotation=15, ha="right")
    ax0.set_ylabel("median relative error (%)")
    ax0.set_title("A. Synthetic stress tests", loc="left")
    ax0.grid(axis="y", alpha=0.22)
    ax0.legend(frameon=False, fontsize=10)
    ax0.text(0.02, 0.94, "2nd order improves\nweak-to-strong medians", transform=ax0.transAxes, fontsize=11, color=BLUE)

    ax1 = fig.add_subplot(gs[0, 1])
    reversals = int(surrogate["braess_like_count"])
    total = int(surrogate["n_lines_tested"])
    ax1.pie(
        [reversals, total - reversals],
        colors=[RED, "#D9E2EC"],
        startangle=90,
        wedgeprops=dict(width=0.36, edgecolor="white", linewidth=2.0),
    )
    ax1.text(0, 0.07, f"{reversals}/{total}", ha="center", va="center", fontsize=32, fontweight="bold", color=RED)
    ax1.text(0, -0.20, "line reinforcements\nraise $\\nu_c$ but lower $\\zeta_{\\min}$", ha="center", va="center", fontsize=12)
    ax1.set_title("B. IEEE39-derived surrogate", loc="left")

    ax2 = fig.add_subplot(gs[0, 2])
    labels = ["diagonal", "second\norder", "reduced\nQEP", "adaptive"]
    med = np.array(
        [
            surrogate["errors"]["diag"]["median"] * 100.0,
            surrogate["errors"]["second"]["median"] * 100.0,
            surrogate["errors"]["rqep"]["median"] * 100.0,
            surrogate["errors"]["adaptive"]["median"] * 100.0,
        ]
    )
    colors = ["#87AFCB", BLUE, GREEN, CYAN]
    ax2.bar(labels, med, color=colors)
    ax2.set_yscale("log")
    ax2.set_ylabel("median error (%)")
    ax2.set_title("C. Surrogate estimator error", loc="left")
    ax2.grid(axis="y", alpha=0.22)
    for i, value in enumerate(med):
        ax2.text(i, value * 1.12, f"{value:.3f}%", ha="center", fontsize=10, fontweight="bold")
    ax2.tick_params(axis="x", labelsize=9)
    ax2.text(0.98, 0.05, "lower is better", transform=ax2.transAxes, ha="right", fontsize=9, color=GRAY)

    fig.suptitle("Evidence is deliberately staged: mechanism first, benchmark claims pending", fontsize=20, fontweight="bold", y=1.02)
    save(fig, "fig_evidence_compact")


def main() -> None:
    setup_style()
    case = load_case()
    plot_witness_bar_ranking(case)
    plot_reversal_story(case)
    plot_method_pipeline()
    plot_evidence_compact(case)
    print("Generated enhanced poster assets:")
    for name in ["fig_witness_bar_ranking", "fig_reversal_story", "fig_method_pipeline", "fig_evidence_compact"]:
        print(FIG / f"{name}.pdf")


if __name__ == "__main__":
    main()
