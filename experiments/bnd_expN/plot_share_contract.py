"""Plot the declared generator-component share contract."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def main() -> None:
    with (OUT / "TABLE_N02_rho_semantics.csv").open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
    for bus, color in ((31, "#2563eb"), (38, "#d97706"), (39, "#059669")):
        subset = [r for r in rows if int(r["bus"]) == bus]
        rho = [float(r["rho"]) for r in subset]
        for ax, quantity, unit in ((axes[0], "P", "MW"), (axes[1], "Q", "Mvar")):
            sg = [float(r[f"SG_target_{quantity}_{unit}"]) for r in subset]
            gf = [float(r[f"GFL_target_{quantity}_{unit}"]) for r in subset]
            ax.plot(rho, sg, "--", color=color, label=f"SG bus {bus}")
            ax.plot(rho, gf, "-", color=color, label=f"GFL bus {bus}")
    axes[0].set_ylabel("Generator active power (MW)")
    axes[1].set_ylabel("Generator reactive power (MVAr)")
    for ax in axes:
        ax.set_xlabel("Converted fraction ρ")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "FIG_N01_PQ_share_vs_rho.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
