"""F7 figures: policy-dependent incompatibility hypergraph maps.

Fill: kappa_Gamma as an ordinal one-hue ramp (validated with the dataviz
validator, --ordinal), no-destabilizing-coalition in neutral grey, base-unstable
hatched. Marks: the exact located crossings, black where kappa changes, white
where H changes at constant kappa, orange where the crossing is a tangency.
Text: the kappa-witness coalitions of the largest H regions.
"""

from __future__ import annotations

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import ndimage  # noqa: E402

from _bootstrap import RESULTS, ROOT  # noqa: E402
from F7_report import BAD, kappa_of, points_file, raster, witnesses  # noqa: E402

F7 = RESULTS / "F7"
FIG = ROOT / "figures"
Y = {
    "F7A": ("k", "excitation-gain scale  k  (t = 1.5)"),
    "F7B": ("t", "excitation time-constant scale  t  (k = 1)"),
    "F7C": ("h", "excitation heterogeneity  h_E  (k = t = 1)"),
}

INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
NONE_GREY = "#f0efec"
UNSTABLE_GREY = "#c3c2b7"
RAMP = {4: "#86b6ef", 3: "#3987e5", 2: "#1c5cab", 1: "#0d366b"}
TANGENCY = "#eb6834"
CODES = {-1: 0, 4: 1, 3: 2, 2: 3, 1: 4, -2: 5}
CMAP = ListedColormap([NONE_GREY, RAMP[4], RAMP[3], RAMP[2], RAMP[1], UNSTABLE_GREY])


def short(label: str) -> str:
    return label.replace("+", "·").replace("|", ", ")


BAND_EDGE = "#e87ba4"


def legend_handles(tongue=False, band_edge=False):
    handles = [
        Patch(fc=NONE_GREY, ec=MUTED, lw=0.3, label="no destabilizing coalition"),
        *[Patch(fc=RAMP[k], label=f"κ = {k}") for k in (4, 3, 2, 1)],
        Patch(
            fc="white", ec=MUTED, hatch="////", lw=0.3, label="base unstable (excluded)"
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            ms=2,
            color=INK,
            label="κ boundary (imaginary axis)",
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            ms=2,
            mfc=SURFACE,
            mec=MUTED,
            mew=0.3,
            label="H boundary at constant κ",
        ),
    ]
    if band_edge:
        handles.append(
            Line2D(
                [],
                [],
                ls="none",
                marker="o",
                ms=4,
                mfc="none",
                mec=BAND_EDGE,
                mew=0.8,
                label="band-edge (0.3 Hz) classification boundary",
            )
        )
    if tongue:
        handles += [
            Line2D(
                [],
                [],
                ls="none",
                marker="v",
                ms=5,
                mfc=TANGENCY,
                mec=INK,
                mew=0.5,
                label="fold: tip of instability tongue",
            ),
            Line2D(
                [],
                [],
                ls="none",
                marker="^",
                ms=5,
                mfc=TANGENCY,
                mec=INK,
                mew=0.5,
                label="fold: closing point of stable gap",
            ),
            Line2D(
                [],
                [],
                ls=(0, (1, 2)),
                lw=0.6,
                color=INK,
                label="pure-policy path of the inset",
            ),
        ]
    return handles


def draw(
    ax,
    name,
    *,
    annotate=8,
    min_fraction=0.004,
    tongue=None,
    legend=True,
    avoid=None,
    short_x=False,
):
    y_name, y_title = Y[name]
    points = pd.read_csv(points_file(name), low_memory=False)
    leaves = pd.read_csv(F7 / f"{name}_leaves.csv")
    events = pd.read_csv(F7 / f"{name}_events.csv", low_memory=False)
    grid = raster(points, leaves)
    xs = (
        points.groupby("I")
        .x.first()
        .reindex(range(grid.shape[0]))
        .interpolate()
        .to_numpy()
    )
    ys = (
        points.groupby("J")
        .y.first()
        .reindex(range(grid.shape[1]))
        .interpolate()
        .to_numpy()
    )
    codes = np.vectorize(lambda s: CODES.get(kappa_of(s), 5) if s else 5)(grid)
    ax.pcolormesh(
        np.sqrt(xs),
        ys,
        codes.T,
        cmap=CMAP,
        vmin=-0.5,
        vmax=5.5,
        shading="nearest",
        rasterized=True,
    )
    unstable = (codes == 5).T
    if unstable.any():
        ax.contourf(
            np.sqrt(xs),
            ys,
            unstable.astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["////"],
        )

    sub = events[(events.kind == "SUBSET_CROSSING") & events.changes_H.astype(bool)]
    gx = np.sqrt(sub.g_star.clip(lower=0))
    gy = sub[f"{y_name}_star"]
    kchange = sub.kappa_before != sub.kappa_after
    tang = sub.boundary_type == "IMAGINARY_AXIS_TANGENCY"
    ax.scatter(
        gx[~kchange & ~tang],
        gy[~kchange & ~tang],
        s=0.15,
        c=SURFACE,
        lw=0,
        rasterized=True,
    )
    ax.scatter(
        gx[kchange & ~tang], gy[kchange & ~tang], s=0.25, c=INK, lw=0, rasterized=True
    )
    ax.scatter(gx[tang], gy[tang], s=6, c=TANGENCY, marker="x", lw=0.6)
    other = sub.boundary_type.astype(str).str.contains("EDGE|CORNER")
    if other.any():
        ax.scatter(
            gx[other], gy[other], s=8, facecolors="none", edgecolors=BAND_EDGE, lw=0.6
        )

    # annotate the largest H regions at the point deepest inside their largest component
    labels = [s for s in set(grid.ravel()) if s and s not in BAD]
    sizes = {s: (grid == s).sum() for s in labels}
    total = sum(sizes.values())
    placed = 0
    for label in sorted(labels, key=lambda s: -sizes[s]):
        if placed >= annotate or sizes[label] / total < min_fraction:
            break
        mask = grid == label
        if avoid is not None:
            gx0, gx1 = np.sqrt(xs[0]), np.sqrt(xs[-1])
            ax_x = (np.sqrt(xs) - gx0) / (gx1 - gx0)
            ax_y = (ys - ys[0]) / (ys[-1] - ys[0])
            blocked = (
                (ax_x[:, None] >= avoid[0])
                & (ax_x[:, None] <= avoid[0] + avoid[2])
                & (ax_y[None, :] >= avoid[1])
                & (ax_y[None, :] <= avoid[1] + avoid[3])
            )
            mask = mask & ~blocked
            if not mask.any():
                continue
        comp, n = ndimage.label(mask)
        biggest = 1 + int(np.argmax(ndimage.sum(mask, comp, range(1, n + 1))))
        # pad with background so labels stay clear of the axes frame
        depth = ndimage.distance_transform_edt(np.pad(comp == biggest, 12))[
            12:-12, 12:-12
        ]
        i, j = np.unravel_index(np.argmax(depth), depth.shape)
        k = kappa_of(label)
        text = (
            "no coalition" if label == "EMPTY" else f"κ={k}: {short(witnesses(label))}"
        )
        dark = k in (2, 1)
        ax.text(
            np.sqrt(xs[i]),
            ys[j],
            text,
            fontsize=6.2,
            ha="center",
            va="center",
            color=SURFACE if dark else INK,
        )
        placed += 1

    if tongue:
        for key, marker in (("tip_fold", "v"), ("gap_closure_fold", "^")):
            f = tongue.get(key)
            if f:
                ax.plot(
                    np.sqrt(f["g"]),
                    f["k"],
                    marker,
                    ms=5,
                    mfc=TANGENCY,
                    mec=INK,
                    mew=0.5,
                )
        ax.axhline(tongue.get("k_ref", 1.30), color=INK, lw=0.5, ls=(0, (1, 2)))
    ticks = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    ax.set_xticks(np.sqrt(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel(
        "reactive-policy gain  g  (√ scale)"
        if short_x
        else "reactive-policy gain  g  (√ scale; 0 = fixed Q, "
        "1 = E30 voltage regulation)"
    )
    ax.set_ylabel(y_title)
    ax.set_title(name, loc="left", fontsize=9, color=MUTED)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color(MUTED)
    ax.tick_params(width=0.5, color=MUTED, labelsize=7)
    if legend:
        ax.legend(
            handles=legend_handles(bool(tongue), bool(other.any())),
            fontsize=6.5,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=False,
        )
    return bool(other.any())


def tongue_inset(ax, tongue):
    grid = pd.read_csv(F7 / "F7_tongue_R_grid.csv", index_col="g")
    col = f"{tongue['k_ref']:.3f}"
    r = grid[col]
    ax.axhline(0, color=MUTED, lw=0.5)
    ax.set_ylim(-0.1, 0.32)
    ax.plot(np.sqrt(r.index), r.values, color=RAMP[2], lw=1.6)
    for c in tongue["crossings_at_k_ref"]:
        ax.plot(np.sqrt(c["g"]), 0, "o", ms=3, color=INK)
    # kappa segments along the same line, read off the map itself
    points = pd.read_csv(points_file("F7A"), usecols=["g", "k", "label"])
    line = points[
        np.isclose(points.k, points.k.iloc[(points.k - tongue["k_ref"]).abs().argmin()])
    ]
    line = line[line.g <= 0.3].sort_values("g")
    kap = [kappa_of(s) for s in line.label]
    start = 0
    for i in range(1, len(kap) + 1):
        if i == len(kap) or kap[i] != kap[start]:
            mid = np.sqrt(0.5 * (line.g.iloc[start] + line.g.iloc[i - 1]))
            text = "∞" if kap[start] == -1 else f"κ={kap[start]}"
            ax.text(
                mid,
                0.27 if kap[start] != -1 else -0.075,
                text,
                fontsize=6,
                ha="center",
                color=INK,
            )
            start = i
    ticks = [0, 0.01, 0.05, 0.1, 0.2]
    ax.set_xticks(np.sqrt(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=6)
    ax.set_xlim(0, np.sqrt(0.3))
    ax.tick_params(labelsize=6, width=0.4)
    ax.set_ylabel("Re λ  (1/s)", fontsize=6)
    ax.set_title(
        f"inter-area family along g at k = {tongue['k_ref']:.2f}",
        fontsize=6.5,
        loc="left",
    )


def main(argv) -> int:
    FIG.mkdir(exist_ok=True)
    tongue = json.loads((F7 / "F7_tongue.json").read_text(encoding="utf-8"))
    names = [n for n in Y if (F7 / f"{n}_events.csv").exists()]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.facecolor": SURFACE,
            "figure.facecolor": "white",
        }
    )
    for name in names:
        fig, ax = plt.subplots(figsize=(8.8, 5.0), dpi=200)
        box = [0.60, 0.06, 0.37, 0.28]
        draw(
            ax,
            name,
            tongue=tongue if name == "F7A" else None,
            avoid=[box[0] - 0.06, 0.0, box[2] + 0.1, box[3] + 0.12]
            if name == "F7A"
            else None,
        )
        if name == "F7A":
            inset = ax.inset_axes(box)
            tongue_inset(inset, tongue)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(FIG / f"{name}_hypergraph_map.{ext}")
        plt.close(fig)
    if len(names) == 3:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), dpi=200)
        edges = False
        for ax, name in zip(axes, names, strict=True):
            edges |= draw(
                ax,
                name,
                annotate=4,
                min_fraction=0.02,
                tongue=tongue if name == "F7A" else None,
                legend=False,
                short_x=True,
            )
        fig.legend(
            handles=legend_handles(True, edges),
            fontsize=7,
            loc="lower center",
            ncol=6,
            frameon=False,
            bbox_to_anchor=(0.5, 0.0),
        )
        fig.tight_layout(rect=(0, 0.1, 1, 1))
        for ext in ("png", "pdf"):
            fig.savefig(FIG / f"F7_three_slices.{ext}")
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
