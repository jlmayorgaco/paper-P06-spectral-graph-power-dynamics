"""F12 figure: Kundur kappa maps, band-limited Gamma and full right half plane."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from _bootstrap import RESULTS, ROOT  # noqa: E402
from F7_figures import (  # noqa: E402
    CMAP,
    MUTED,
    NONE_GREY,
    RAMP,
    UNSTABLE_GREY,
)
from F7_report import kappa_of  # noqa: E402

CODES = {-1: 0, 4: 1, 3: 2, 2: 3, 1: 4, -2: 5}


def panel(ax, frame, column, y, title):
    grid = frame.pivot_table(index="j", columns="i", values=column, aggfunc="first")
    codes = np.vectorize(
        lambda s: CODES.get(
            kappa_of(s) if s not in ("BASE_UNSTABLE", "INFEASIBLE") else -2, 5
        )
    )(grid.to_numpy())
    g = frame.groupby("i").g.first().to_numpy()
    yy = frame.groupby("j")[y].first().to_numpy()
    ax.pcolormesh(
        np.sqrt(g),
        yy,
        codes,
        cmap=CMAP,
        vmin=-0.5,
        vmax=5.5,
        shading="nearest",
        rasterized=True,
    )
    ticks = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    ax.set_xticks(np.sqrt(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel("reactive-policy gain  g  (sqrt scale)")
    ax.set_ylabel({"k": "SEXS gain scale  k", "t": "SEXS time-constant scale  t"}[y])
    ax.set_title(title, loc="left", fontsize=9, color=MUTED)
    ax.tick_params(labelsize=7)


def main() -> int:
    nodes = pd.read_csv(RESULTS / "F12" / "F12_nodes.csv")
    full = pd.read_csv(RESULTS / "F12" / "F12_full_rhp.csv")
    nodes = nodes.merge(
        full[["map", "i", "j", "label_full"]], on=["map", "i", "j"], how="left"
    )
    nodes["label_full"] = nodes.label_full.fillna(nodes.label)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=180)
    for row, (m, y) in enumerate((("K12A", "k"), ("K12B", "t"))):
        f = nodes[nodes["map"] == m]
        panel(
            axes[row, 0],
            f,
            "label",
            y,
            f"{m}: protected band 0.3-1.5 Hz (preregistered)",
        )
        panel(
            axes[row, 1], f, "label_full", y, f"{m}: whole right half plane (post-hoc)"
        )
    handles = [
        Patch(fc=NONE_GREY, ec=MUTED, lw=0.3, label="no destabilizing coalition"),
        *[Patch(fc=RAMP[k], label=f"kappa = {k}") for k in (4, 3, 2, 1)],
        Patch(fc=UNSTABLE_GREY, label="base unstable"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"F12_kundur_maps.{ext}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
