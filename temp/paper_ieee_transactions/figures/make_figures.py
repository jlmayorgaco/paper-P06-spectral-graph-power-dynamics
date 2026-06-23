from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def damping_region():
    rng = np.random.default_rng(12)
    nu = np.linspace(0.2, 12, 300)
    zeta_star = 0.05
    boundary = 2 * zeta_star * np.sqrt(nu)

    fig, axes = plt.subplots(1, 2, figsize=(7.05, 2.35), sharey=True)

    pts_nu = np.exp(rng.uniform(np.log(0.25), np.log(12), 34))
    pts_delta = 0.07 + 0.22 * rng.random(34) * np.sqrt(pts_nu) / np.sqrt(12)
    unsafe = pts_delta < 2 * zeta_star * np.sqrt(pts_nu)

    ax = axes[0]
    ax.plot(nu, boundary, color="#111111", lw=1.8)
    ax.fill_between(nu, boundary, boundary.max() * 1.35, color="#d9ead3", alpha=0.75)
    ax.fill_between(nu, 0, boundary, color="#f4cccc", alpha=0.75)
    ax.scatter(pts_nu[~unsafe], pts_delta[~unsafe], s=16, color="#2b7a0b", label="safe")
    ax.scatter(pts_nu[unsafe], pts_delta[unsafe], s=24, color="#b00020", marker="x", label="unsafe")
    ax.set_title("commuting damping", fontsize=8)
    ax.set_xlabel(r"modal stiffness $\nu_k$")
    ax.set_ylabel(r"modal damping $\delta_k$")
    ax.legend(frameon=False, fontsize=7, loc="upper left")

    ax = axes[1]
    ax.plot(nu, boundary, color="#111111", lw=1.8, label=r"$2\zeta_\star\sqrt{\nu}$")
    ax.fill_between(nu, boundary, boundary.max() * 1.35, color="#d9ead3", alpha=0.75)
    ax.fill_between(nu, 0, boundary, color="#f4cccc", alpha=0.75)
    # A stylized coupled case: the diagonal point is above the target boundary,
    # but the QEP-corrected margin moves below the boundary.
    nu0 = 7.6
    delta0 = 2 * zeta_star * np.sqrt(nu0) + 0.035
    nu1 = 7.6
    delta1 = 2 * zeta_star * np.sqrt(nu1) - 0.030
    ax.scatter([nu0], [delta0], s=32, color="#2b7a0b", label="diagonal screen")
    ax.scatter([nu1], [delta1], s=42, color="#b00020", marker="x", label="coupled QEP")
    ax.annotate("", xy=(nu1, delta1), xytext=(nu0, delta0), arrowprops=dict(arrowstyle="->", lw=1.0))
    bg_nu = np.exp(rng.uniform(np.log(0.25), np.log(12), 22))
    bg_delta = 0.08 + 0.18 * rng.random(22) * np.sqrt(bg_nu) / np.sqrt(12)
    ax.scatter(bg_nu, bg_delta, s=12, color="#8c8c8c", alpha=0.45)
    ax.set_title("non-proportional coupling", fontsize=8)
    ax.set_xlabel(r"modal stiffness $\nu_k$")
    ax.legend(frameon=False, fontsize=7, loc="upper left")

    for ax in axes:
        ax.set_xlim(0, 12)
        ax.set_ylim(0, boundary.max() * 1.25)
        ax.grid(True, alpha=0.22)
    save(fig, "fig1_damping_region")


def estimator_errors():
    labels = ["prop.", "weak", "moderate", "strong", "clustered"]
    diag = np.array([0.0, 0.03, 2.23, 5.05, 300.2])
    qep = np.array([0.0, 0.018, 1.31, 4.00, 149.7])
    x = np.arange(len(labels))
    w = 0.36
    fig, ax = plt.subplots(figsize=(3.45, 2.35))
    ax.bar(x - w / 2, diag + 1e-3, width=w, color="#7aa6c2", label="diagonal")
    ax.bar(x + w / 2, qep + 1e-3, width=w, color="#d08c60", label="reduced QEP")
    ax.set_yscale("log")
    ax.set_ylabel("median relative error (%)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=7)
    ax.text(0.02, 0.98, "synthetic stress tests", transform=ax.transAxes, va="top", fontsize=7)
    save(fig, "fig2_estimator_errors")


def workflow():
    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    ax.axis("off")
    boxes = [
        (0.5, 0.88, "Build\n$\\widetilde L,\\widetilde D$"),
        (0.5, 0.68, "Compute\n$\\nu_k,q_k,\\Gamma$"),
        (0.5, 0.48, "Diagonal screen\n$\\min_k\\delta_k/(2\\sqrt{\\nu_k})$"),
        (0.22, 0.25, "weak coupling:\naccept diagonal"),
        (0.78, 0.25, "cluster/coupling:\nreduced QEP"),
        (0.5, 0.06, "planning scores:\nconversion, damping, lines"),
    ]
    for x, y, text in boxes:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.28", fc="#f8f8f8", ec="#444444", lw=0.8),
        )
    arrow = dict(arrowstyle="->", lw=0.9, color="#333333")
    for y1, y2 in [(0.82, 0.74), (0.62, 0.54)]:
        ax.annotate("", xy=(0.5, y2), xytext=(0.5, y1), arrowprops=arrow)
    ax.annotate("", xy=(0.22, 0.32), xytext=(0.43, 0.42), arrowprops=arrow)
    ax.annotate("", xy=(0.78, 0.32), xytext=(0.57, 0.42), arrowprops=arrow)
    ax.annotate("", xy=(0.43, 0.11), xytext=(0.25, 0.19), arrowprops=arrow)
    ax.annotate("", xy=(0.57, 0.11), xytext=(0.75, 0.19), arrowprops=arrow)
    save(fig, "fig3_adaptive_workflow")


def line_tradeoff():
    omega = np.linspace(0.5, 5.5, 300)
    sigma_const = -0.12
    zeta_const = -sigma_const / np.sqrt(sigma_const**2 + omega**2)
    sigma_improved = -(0.12 + 0.025 * omega)
    zeta_improved = -sigma_improved / np.sqrt(sigma_improved**2 + omega**2)

    fig, ax = plt.subplots(figsize=(3.45, 2.35))
    ax.plot(omega, zeta_const, color="#b00020", lw=1.8, label="frequency rises, decay fixed")
    ax.plot(omega, zeta_improved, color="#2b7a0b", lw=1.8, label="decay rises enough")
    ax.set_xlabel(r"oscillation frequency $|\omega|$")
    ax.set_ylabel(r"damping ratio $\zeta$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=7)
    ax.text(0.04, 0.95, r"$\zeta=-\sigma/\sqrt{\sigma^2+\omega^2}$", transform=ax.transAxes, va="top", fontsize=8)
    save(fig, "fig4_line_tradeoff")


if __name__ == "__main__":
    damping_region()
    estimator_errors()
    workflow()
    line_tradeoff()
