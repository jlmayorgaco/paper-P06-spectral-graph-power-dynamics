"""Render the post-freeze PowerDynamics rho/PLL validation grid."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "reports" / "experiment_D" / "tables" / "TABLE_D17_PD_rho_PLL_grid.csv"
OUT = ROOT / "reports" / "experiment_D" / "figures"
SIGMA_REQ = 0.05


def main() -> None:
    data = pd.read_csv(TABLE)
    expected_buses = [30, 33, 35, 37]
    if len(data) != 132 or sorted(data["bus"].unique().tolist()) != expected_buses:
        raise RuntimeError("Expected the complete 4-bus, 11-rho, 3-beta validation grid")

    data["sigma_margin"] = -data["alpha"] - SIGMA_REQ
    vmin = min(-0.45, float(data["sigma_margin"].min()))
    vmax = max(0.10, float(data["sigma_margin"].max()))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.7), constrained_layout=True, sharex=True, sharey=True)
    image = None
    for ax, bus in zip(axes.flat, expected_buses):
        subset = data.loc[data["bus"] == bus]
        matrix = subset.pivot(index="beta", columns="rho", values="sigma_margin")
        x = matrix.columns.to_numpy(dtype=float)
        y = matrix.index.to_numpy(dtype=float)
        image = ax.imshow(
            matrix.to_numpy(dtype=float),
            origin="lower",
            interpolation="nearest",
            aspect="auto",
            extent=(x.min() - 0.05, x.max() + 0.05, y.min() - 0.05, y.max() + 0.05),
            cmap="RdYlGn",
            norm=norm,
        )
        feasible_count = int(subset["within_sigma"].astype(bool).sum())
        ax.set_title(f"Bus {bus}: {feasible_count}/33 grid cells meet the margin")
        ax.set_xticks(x)
        ax.set_yticks(y)
        ax.grid(which="major", color="white", linewidth=0.7, alpha=0.75)
        for _, row in subset.loc[~subset["within_sigma"].astype(bool)].iterrows():
            ax.text(row["rho"], row["beta"], "×", ha="center", va="center", color="black", fontsize=14, weight="bold")
        ax.set_xlabel(r"Replacement share $\rho$")
        ax.set_ylabel(r"PLL scale $\beta$")

    colorbar = fig.colorbar(image, ax=axes, shrink=0.91, pad=0.025)
    colorbar.set_label(r"Spectral margin $-\alpha-0.05$ (s$^{-1}$); green = feasible")
    fig.suptitle("Post-freeze PowerDynamics validation grid", fontsize=14, weight="bold")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "FIG_D01_PD_validation_grid.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / "FIG_D01_PD_validation_grid.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
