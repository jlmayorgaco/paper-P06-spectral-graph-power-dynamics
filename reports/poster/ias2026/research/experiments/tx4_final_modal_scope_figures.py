"""Figures for the TX4 final modal-scope correction audit."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
OUT = ROOT / "figures" / "final_audit"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.dpi": 140,
        "savefig.dpi": 260,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=260, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def local_vs_collective(local: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    colors = {30: "#277da1", 33: "#f9844a", 35: "#90be6d", 37: "#f9c74f"}
    for bus, frame in local.groupby("device_i"):
        frame = frame.sort_values("g")
        ax.semilogy(frame.g, frame.local_sigma_min_I_plus_Mii, "o-", color=colors[int(bus)], label=f"local I+M_{{{int(bus)}{int(bus)}}}")
    unique = local.drop_duplicates("g").sort_values("g")
    ax.semilogy(unique.g, unique.collective_sigma_min_I_plus_QH, "k--", lw=2.2, label="collective I+Q_H")
    g_star = 0.20768140519037842
    ax.axvline(g_star, color="#d95f02", ls=":", lw=1.2, label="H4 boundary")
    ax.set_xlabel("controller gain g")
    ax.set_ylabel("smallest singular value")
    ax.set_title("Physical local factors stay regular while collective closure closes")
    ax.set_ylim(1e-9, 2)
    ax.legend(fontsize=8, ncol=2)
    save(fig, "F_TRUE_LOCAL_VS_COLLECTIVE")


def sign_and_return(local: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
    q = local.drop_duplicates("g").sort_values("g")
    g_star = 0.20768140519037842
    axes[0].plot(q.q_eigenvalue_real, q.q_eigenvalue_imag, "o-", color="#277da1")
    axes[0].scatter([-1], [0], marker="x", s=70, color="#d95f02", label="target -1")
    axes[0].set_title("Q_H eigenvalue")
    axes[0].set_xlabel("real")
    axes[0].set_ylabel("imaginary")
    axes[0].legend(fontsize=8)
    for bus, frame in local.groupby("device_i"):
        frame = frame.sort_values("g")
        axes[1].plot(frame.return_eigenvalue_real, frame.return_eigenvalue_imag, "o-", label=f"bus {int(bus)}")
    axes[1].scatter([1], [0], marker="x", s=70, color="#d95f02", label="target +1")
    axes[1].set_title("contextual return eigenvalue")
    axes[1].set_xlabel("real")
    axes[1].set_ylabel("imaginary")
    axes[1].legend(fontsize=7)
    axes[2].plot(q.g, q.alpha_EM, "o-", color="#277da1", label="alpha_EM")
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].axvline(g_star, color="#d95f02", ls=":", label="boundary")
    axes[2].set_title("critical EM-band alpha")
    axes[2].set_xlabel("controller gain g")
    axes[2].set_ylabel("alpha [s^-1]")
    axes[2].legend(fontsize=8)
    fig.suptitle("Correct signs: Q_H reaches -1; return reaches +1")
    save(fig, "F_Q_MINUS1_RETURN_PLUS1")


def tds_equal_horizon(root: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    specs = [
        ("trace_proper_30p33p35_D2.csv.gz", "proper 30+33+35", "#1b9e77"),
        ("trace_H4_original_30p33p35p37_D2.csv.gz", "H4 original", "#d95f02"),
        ("trace_H4_retuned_g_0.25_D2.csv.gz", "H4 retuned g=0.25", "#277da1"),
    ]
    for filename, label, color in specs:
        path = root / "results" / "TX4_TDS" / filename
        trace = pd.read_csv(path)
        trace = trace[trace.t <= 30.0]
        envelope = np.abs(trace.signal - trace.signal.mean())
        ax.plot(trace.t, np.maximum(envelope, 1e-10), color=color, lw=1.1, label=label)
    ax.axvspan(1.0, 1.1, color="#888888", alpha=0.18, label="common pulse")
    ax.set_yscale("log")
    ax.set_xlim(0, 30)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("absolute speed-difference deviation")
    ax.set_title("Equal-horizon TDS comparison: common +2% bus-20 pulse")
    ax.legend(fontsize=8)
    save(fig, "F_TDS_FINAL_EQUAL_HORIZON")


def v4_lattice(v4: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    levels = {k: v4[v4.cardinality == k].sort_values("portfolio").reset_index(drop=True) for k in range(5)}
    pos: dict[str, tuple[float, float]] = {}
    for k, frame in levels.items():
        n = len(frame)
        for j, row in frame.iterrows():
            x = (j - (n - 1) / 2) * (0.86 if k in (1, 3) else 0.72) + 2.0
            y = float(k)
            pos[row.portfolio] = (x, y)
    for k in range(4):
        for child in levels[k + 1].portfolio:
            members = set(child.split("+"))
            for parent in levels[k].portfolio:
                if set(parent.split("+")) <= members and len(set(parent.split("+"))) == k:
                    ax.plot([pos[parent][0], pos[child][0]], [pos[parent][1], pos[child][1]], color="#bdbdbd", lw=0.7, zorder=1)
    for portfolio, (x, y) in pos.items():
        row = v4[v4.portfolio == portfolio].iloc[0]
        unstable = row.truth_verdict == "UNSTABLE_PREDICTED"
        color = "#d95f02" if unstable else "#1b9e77"
        ax.scatter([x], [y], s=170, color=color, edgecolor="black", zorder=3)
        if portfolio == "30+33+35+37":
            ax.text(x, y, "H4", ha="center", va="center", color="white", weight="bold", fontsize=8)
        elif portfolio == "BASE":
            ax.text(x, y, "all SG", ha="center", va="center", color="white", weight="bold", fontsize=7)
    h4 = v4[v4.portfolio == "30+33+35+37"].iloc[0]
    proper = v4[v4.portfolio != "30+33+35+37"]
    ax.annotate(f"H4: alpha={h4.alpha_truth:+.3f} s^-1, f={h4.frequency_truth_hz:.3f} Hz", xy=pos["30+33+35+37"], xytext=(4.8, 3.5), arrowprops=dict(arrowstyle="->", color="#d95f02"), color="#d95f02", fontsize=8)
    ax.text(4.8, 2.85, f"worst proper alpha={proper.alpha_truth.max():+.3f} s^-1", color="#1b9e77", fontsize=8)
    ax.set_xlim(-2.0, 6.0)
    ax.set_ylim(-0.45, 4.55)
    ax.set_yticks(range(5), ["0", "1", "2", "3", "4"])
    ax.set_xlabel("replacement cardinality")
    ax.set_ylabel("lattice level")
    ax.set_title("V4 full-order lattice: H4 is the sole true blocker")
    save(fig, "F_V4_TRUE_LATTICE")


def main() -> None:
    root = ROOT
    local = pd.read_csv(root / "results" / "TX4_TRUE_LOCAL_VS_COLLECTIVE.csv")
    v4 = pd.read_csv(root / "results" / "TX4_V4_BLIND_VS_FULL.csv")
    local_vs_collective(local)
    sign_and_return(local)
    tds_equal_horizon(root)
    v4_lattice(v4)
    print(f"created figures in {OUT}")


if __name__ == "__main__":
    main()
