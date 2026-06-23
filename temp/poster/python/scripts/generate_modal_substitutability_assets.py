from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
POSTER = ROOT / "poster_ieee_session"
ASSET_DIR = POSTER / "modal_substitutability_assets"
FIG = POSTER / "figures"

DARK = "#16251d"
GREEN = "#18884a"
RED = "#bd3a31"
GRAY = "#8d9690"
LIGHT = "#f7fbf8"


def load_data() -> dict:
    with (ASSET_DIR / "poster_data.json").open("r", encoding="utf-8") as f:
        return json.load(f)


def save(fig: plt.Figure, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=360, bbox_inches="tight")
    plt.close(fig)


def build_laplacian(edges: list[list[float]], n: int) -> np.ndarray:
    lap = np.zeros((n, n), dtype=float)
    for i, j, weight in edges:
        i = int(i)
        j = int(j)
        weight = float(weight)
        lap[i, j] -= weight
        lap[j, i] -= weight
        lap[i, i] += weight
        lap[j, j] += weight
    return lap


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": DARK,
            "axes.labelcolor": DARK,
            "xtick.color": DARK,
            "ytick.color": DARK,
            "text.color": DARK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def plot_substitutability(data: dict) -> None:
    subst = np.array(data["substitutability"], dtype=float)
    edges = data["edges"]
    best = int(data["best_node"])
    worst = int(data["worst_node"])

    pos = {
        0: (0.00, 1.08),
        1: (1.03, 1.78),
        2: (1.02, 0.34),
        3: (2.13, 1.78),
        4: (2.12, 0.34),
        5: (3.18, 1.08),
    }

    fig, ax = plt.subplots(figsize=(12.4, 8.1))
    ax.set_facecolor(LIGHT)

    max_w = max(float(edge[2]) for edge in edges)
    for i, j, weight in edges:
        x = [pos[int(i)][0], pos[int(j)][0]]
        y = [pos[int(i)][1], pos[int(j)][1]]
        ax.plot(x, y, color=GRAY, lw=2.0 + 3.2 * float(weight) / max_w, alpha=0.62, zorder=1)

    cmap = plt.cm.RdYlGn
    for node, value in enumerate(subst):
        edge_color = DARK
        line_width = 3.2
        size = 3800
        if node in {best, worst}:
            line_width = 5.0
            size = 4300
        ax.scatter(
            *pos[node],
            s=size,
            color=cmap(value),
            edgecolors=edge_color,
            linewidths=line_width,
            zorder=3,
        )
        ax.text(
            *pos[node],
            f"Node {node}\ns={value:.2f}",
            ha="center",
            va="center",
            fontsize=21,
            fontweight="bold",
            color="black",
            zorder=4,
        )

    ax.text(
        pos[best][0],
        pos[best][1] - 0.34,
        "topology-cheap",
        ha="center",
        va="top",
        fontsize=18,
        fontweight="bold",
        color=GREEN,
    )
    ax.text(
        pos[worst][0],
        pos[worst][1] + 0.34,
        "inertia/SynCon",
        ha="center",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        color=RED,
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("modal substitutability $s_i$", fontsize=20, labelpad=12)
    cbar.ax.tick_params(labelsize=17)

    ax.set_xlim(-0.30, 3.48)
    ax.set_ylim(-0.15, 2.18)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "Six-bus modal substitutability map",
        fontsize=30,
        fontweight="bold",
        pad=18,
    )
    ax.text(
        0.5,
        -0.045,
        "green = reinforce incident lines | red = virtual inertia / synchronous condenser",
        transform=ax.transAxes,
        ha="center",
        fontsize=18,
        color=DARK,
    )
    save(fig, "fig_modal_substitutability")


def plot_opposite_levers(data: dict) -> None:
    inertia = np.array(data["inertia"], dtype=float)
    edges = data["edges"]
    n = len(inertia)
    lap = build_laplacian(edges, n)
    minv = np.diag(1.0 / np.sqrt(inertia))
    ltilde = minv @ lap @ minv
    nu, qmat = np.linalg.eigh(ltilde)
    crit = int(np.argmax(nu))
    q = qmat[:, crit]

    inertia_shift = []
    for i in range(n):
        eye_i = np.zeros((n, n))
        eye_i[i, i] = 1.0
        dlt = -0.5 * (eye_i @ ltilde + ltilde @ eye_i)
        inertia_shift.append(float(q @ dlt @ q))

    line_shift = []
    labels = []
    for i, j, _ in edges:
        b = np.zeros(n)
        b[int(i)] = 1.0
        b[int(j)] = -1.0
        dlt = minv @ np.outer(b, b) @ minv
        line_shift.append(float(q @ dlt @ q))
        labels.append(f"{int(i)}-{int(j)}")

    x_inertia = np.arange(n)
    gap = 1.4
    x_lines = np.arange(len(line_shift)) + n + gap

    fig, ax = plt.subplots(figsize=(13.0, 5.8))
    ax.axhline(0.0, color=DARK, lw=1.8)
    ax.bar(x_inertia, inertia_shift, color=RED, label="virtual inertia at node", width=0.72)
    ax.bar(x_lines, line_shift, color="#2aae62", label="line reinforcement", width=0.72)

    labels_all = [f"n{i}" for i in range(n)] + labels
    ax.set_xticks(list(x_inertia) + list(x_lines))
    ax.set_xticklabels(labels_all, rotation=38, ha="right", fontsize=16)
    ax.tick_params(axis="y", labelsize=17)
    ax.set_ylabel(r"first-order shift of $\nu_c$", fontsize=20)
    ax.set_title("Two levers, opposite signs on the same critical mode", fontsize=27, fontweight="bold", pad=16)
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, fontsize=18, loc="lower right")
    ax.text(
        np.mean(x_inertia),
        min(inertia_shift) * 0.72,
        "inertia lowers $\\nu_c$",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", fc=RED, ec=RED),
    )
    ax.text(
        np.mean(x_lines),
        max(line_shift) * 0.72,
        "lines raise $\\nu_c$",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", fc="#2aae62", ec="#2aae62"),
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    save(fig, "fig_opposite_levers")


def main() -> None:
    setup_style()
    data = load_data()
    plot_substitutability(data)
    plot_opposite_levers(data)
    print("Generated modal substitutability poster assets:")
    print(FIG / "fig_modal_substitutability.pdf")
    print(FIG / "fig_opposite_levers.pdf")


if __name__ == "__main__":
    main()
