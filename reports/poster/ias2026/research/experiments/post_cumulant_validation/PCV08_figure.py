# ruff: noqa: E501  -- plot labels kept on one line
"""PCV08 - Phase 10 four-panel decision figure (frozen data only, no computation).

A  subset lattice of H4 at P4 (PCV02 B0 truth)
B  lower-order hierarchy vs exact verdict (PCV03 hierarchy)
C  clean g-only counterfactual on k = 1.425 (frozen FC03 path) with g* marked
D  boundary motion: exact alpha(H4) along g and k, port tangents, port-guided Newton
   (frozen FC03 path, frozen F7A k-line, PCV04 Newton)
No connected-cumulant panel.
Writes figures/20260911_portfolio_decision_case.{pdf,png,svg} with fixed metadata.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
PCV = RESEARCH / "results" / "PCV"
FC_RUN = RESEARCH / "outputs/ias2026/final_math_nonlinear_validation_20260910T231539"
FIG = RESEARCH / "figures" / "20260911_portfolio_decision_case"
CORE = (30, 33, 35, 37)
G_STAR = 0.20768140450381395
INK, MUTED, GRID = "#1f2328", "#5b6470", "#d0d7de"
UNSTABLE, STABLE = "#b42318", "#1f6f8b"  # status: always with text / marker shape
SERIES = {"g": "#1f6f8b", "k": "#c26a00"}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7.5,
        "axes.titlesize": 8,
        "axes.labelsize": 7.5,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "svg.hashsalt": "pcv08",
        "pdf.fonttype": 42,
    }
)


def label(s):
    return "+".join(map(str, s)) or "BASE"


def panel_a(ax):
    truth = pd.read_csv(PCV / "PCV02" / "PCV02_truth_core.csv")
    t = truth[truth.point == "P4"].set_index("subset")
    pos = {}
    for size in range(5):
        subs = [tuple(s) for s in combinations(CORE, size)]
        half = 0.13 * (len(subs) - 1)
        xs = np.linspace(-half, half, len(subs))
        for x, s in zip(xs, subs, strict=True):
            pos[s] = (x, size)
    for s, (x, y) in pos.items():
        for b in CORE:
            if b not in s:
                sup = tuple(sorted((*s, b)))
                ax.plot(
                    [x, pos[sup][0]], [y, pos[sup][1]], color=GRID, lw=0.6, zorder=1
                )
    for s, (x, y) in pos.items():
        r = t.loc[label(s)]
        unstable = r.status == "UNSTABLE"
        ax.scatter(
            x,
            y,
            s=150 if unstable else 95,
            marker="s" if unstable else "o",
            color=UNSTABLE if unstable else STABLE,
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
        )
        name = "∅" if not s else label(s).replace("+", "·")
        ax.text(x, y + 0.2, name, ha="center", va="bottom", fontsize=5.0, color=INK)
        ax.text(
            x,
            y - (0.3 if unstable else 0.22),
            f"{r.alpha:+.3f}",
            ha="center",
            va="top",
            fontsize=5.0,
            color=UNSTABLE if unstable else MUTED,
        )
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.7, 5.3)
    ax.set_yticks(range(5), ["0", "1", "2", "3", "4"])
    ax.set_ylabel("replaced units |S|")
    ax.set_xticks([])
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.set_title(
        "A  P4 lattice: 15 proper subsets stable, H4 not",
        loc="left",
        color=INK,
    )
    ax.text(
        1.0,
        1.0,
        "α⊥ [s⁻¹] under each node\n■ unstable   ● stable",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color=MUTED,
    )


def panel_b(ax):
    h = pd.read_csv(PCV / "PCV03" / "PCV03_hierarchy.csv")
    p = h[h.point == "P4"]
    exact = float(p[p.level == "exact critical eigenvalue"].predicted_alpha.iloc[0])
    rows = [
        ("α space", "B3 modal sensitivity", "B3 first-order modal sensitivity"),
        ("α space", "B4 additive", "B4 additive singles"),
        ("α space", "B5 pairwise", "B5 pairwise"),
        ("α space", "B5b third order", "B5b third order"),
        ("det space", "det: local only", "local factors only (order <= 1)"),
        ("det space", "det: order ≤ 2", "pairwise closure (order <= 2)"),
        ("det space", "det: order ≤ 3", "third-order closure (order <= 3)"),
        ("exact", "exact closure = B0", "exact full closure (order 4)"),
    ]
    ys = np.arange(len(rows))[::-1]
    xmin = -0.26
    for y, (_space, _name, key) in zip(ys, rows, strict=True):
        r = p[p.level == key].iloc[0]
        v = r.predicted_alpha
        if not np.isfinite(v):
            ax.text(
                xmin + 0.005,
                y,
                "no zero within radius 1",
                va="center",
                fontsize=6,
                color=MUTED,
            )
        else:
            unstable = v > 0
            ax.plot([xmin, v], [y, y], color=GRID, lw=0.8, zorder=1)
            ax.scatter(
                v,
                y,
                s=34,
                marker="s" if unstable else "o",
                color=UNSTABLE if unstable else STABLE,
                edgecolor="white",
                linewidth=0.8,
                zorder=3,
            )
            ax.text(
                v + (0.012 if v < 0.1 else -0.012),
                y + 0.33,
                f"{v:+.3f}",
                ha="left" if v < 0.1 else "right",
                fontsize=5.8,
                color=INK,
            )
    ax.axvline(0, color=INK, lw=0.7)
    ax.axvline(exact, color=UNSTABLE, lw=0.6, ls=":")
    ax.axhline(3.5, color=GRID, lw=0.6)
    ax.axhline(0.5, color=GRID, lw=0.6)
    ax.set_xlim(xmin, 0.19)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(ys, [r[1] for r in rows], fontsize=6.2, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("predicted α⊥(H4) at P4 [s⁻¹]  (> 0 flags H4)")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title(
        "B  Lower-order screens vs exact closure",
        loc="left",
        color=INK,
    )


def panel_c(ax):
    f = pd.read_csv(FC_RUN / "FC03_governed_replication" / "FC03_points.csv")
    path = f[f.tag == "PATH k=1.425"].sort_values("g")
    regimes = [
        (0.0, 0.0128, "H = {30·33·35}\nκ = 3"),
        (0.0128, G_STAR, "H = {H4}, κ = 4"),
        (G_STAR, 1.0, "H = ∅ (all 16 stable)"),
    ]
    shades = ["#f3e8e6", "#fbeeee", "#eef4f7"]
    for (a, b, _text), c in zip(regimes, shades, strict=True):
        ax.axvspan(a, b, color=c, lw=0, zorder=0)
    ax.text(0.11, -0.175, regimes[1][2], ha="center", fontsize=6.0, color=UNSTABLE)
    ax.text(0.60, -0.175, regimes[2][2], ha="center", fontsize=6.0, color=STABLE)
    ax.annotate(
        "κ = 3 for g < 0.013",
        (0.006, 0.305),
        textcoords="offset points",
        xytext=(10, 2),
        fontsize=5.4,
        color=MUTED,
        arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.5},
    )
    ax.plot(path.g, path.alpha_flag_frozen, color=INK, lw=1.4, zorder=2)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.axvline(G_STAR, color=INK, lw=0.7, ls="--")
    ax.text(G_STAR + 0.01, 0.12, f"g* = {G_STAR:.5f}", fontsize=6.2, color=INK)
    cf = pd.read_csv(PCV / "PCV03" / "PCV03_counterfactual.csv").set_index("quantity")
    for pid, g in (("P4", 0.03625), ("G_S", 0.25), ("G_S2", 1.0)):
        a = float(cf.loc["alpha", pid])
        unstable = a > 0
        ax.scatter(
            g,
            a,
            s=36,
            marker="s" if unstable else "o",
            color=UNSTABLE if unstable else STABLE,
            edgecolor="white",
            linewidth=0.9,
            zorder=4,
        )
        ax.annotate(
            f"{pid}\n{a:+.3f}",
            (g, a),
            textcoords="offset points",
            xytext=(6, 6 if unstable else -16),
            fontsize=6,
            color=INK,
        )
    gscr = float(cf.loc["B1/B2/B7 gscr", "P4"])
    pg = float(cf.loc["B1/B2/B7 pg_mw", "P4"])
    sn = float(cf.loc["B1/B2/B7 sn_mva", "P4"])
    ax.text(
        0.99,
        0.97,
        f"identical at every g:\nnetwork, locations, ratings,\nscheduled P/Q, k, t, h,\ngSCR(H4) = {gscr:.3f},\nPg = {pg:.0f} MW, Sn = {sn:.0f} MVA",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color=MUTED,
    )
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.2, 0.34)
    ax.set_xlabel("converter Q/V gain g   (k = 1.425, t = 1.5, h = 1 fixed)")
    ax.set_ylabel("α⊥(H4) [s⁻¹]")
    ax.grid(axis="y", color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title(
        "C  Clean g-only counterfactual (k = 1.425)",
        loc="left",
        color=INK,
    )


def panel_d(ax):
    f = pd.read_csv(FC_RUN / "FC03_governed_replication" / "FC03_points.csv")
    path = f[f.tag == "PATH k=1.425"].sort_values("g")
    a7 = pd.read_csv(RESEARCH / "results/F7/F7A_points.csv.gz", low_memory=False)
    kline = a7[
        np.isclose(a7.g, 0.03625) & (a7.k >= 1.25) & (a7.k <= 1.425)
    ].sort_values("k")
    summ = json.loads(
        (PCV / "PCV04" / "PCV04_summary.json").read_text(encoding="utf-8")
    )
    pol = pd.read_csv(PCV / "PCV04" / "PCV04_policy_task_e.csv").set_index("coordinate")
    a0 = summ["alpha_p4"]
    scale = {"g": 1.0, "k": 1.8}
    x_g = (path.g - 0.03625) / scale["g"]
    x_k = (kline.k - 1.425) / scale["k"]
    ax.plot(
        x_g,
        path.alpha_flag_frozen,
        color=SERIES["g"],
        lw=1.4,
        label="g: exact (frozen FC03 path)",
    )
    ax.plot(
        x_k,
        kline["bandmax_30+33+35+37"],
        color=SERIES["k"],
        lw=0,
        marker="o",
        ms=2.0,
        label="k: exact (frozen F7A grid)",
    )
    spans = {"g": (-0.03, 0.14), "k": (-0.11, 0.03)}
    for c in ("g", "k"):
        xs = np.linspace(*spans[c], 50)
        slope = pol.loc[c, "port_dRe_s"] * scale[c]
        ax.plot(
            xs,
            a0 + slope * xs,
            color=SERIES[c],
            lw=0.8,
            ls=":",
            label=f"{c}: port tangent ({pol.loc[c, 'port_dRe_s']:+.3f})",
        )
        it = summ[f"newton_{c}"]
        xi = [(r[c] - (0.03625 if c == "g" else 1.425)) / scale[c] for r in it]
        ax.scatter(
            xi,
            [r["re_s"] for r in it],
            s=16,
            marker="D",
            color=SERIES[c],
            edgecolor="white",
            linewidth=0.6,
            zorder=4,
            label=f"{c}: port-guided Newton iterates",
        )
    bis = summ["bisect"]
    for c in ("g", "k"):
        xb = (bis[c] - (0.03625 if c == "g" else 1.425)) / scale[c]
        ax.axvline(xb, color=SERIES[c], lw=0.6, ls="--")
    ax.text(
        (bis["g"] - 0.03625) + 0.004,
        -0.17,
        f"g* = {bis['g']:.5f}\nNewton: {len(summ['newton_g']) - 1} steps",
        fontsize=5.8,
        color=SERIES["g"],
    )
    ax.text(
        (bis["k"] - 1.425) / 1.8 - 0.004,
        -0.17,
        f"k* = {bis['k']:.5f}\nNewton: {len(summ['newton_k']) - 1} steps",
        fontsize=5.8,
        color=SERIES["k"],
        ha="right",
    )
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.scatter([0], [a0], s=40, marker="s", color=UNSTABLE, edgecolor="white", zorder=5)
    ax.annotate(
        "P4", (0, a0), textcoords="offset points", xytext=(4, 5), fontsize=6, color=INK
    )
    ax.set_xlim(-0.12, 0.25)
    ax.set_ylim(-0.2, 0.34)
    ax.set_xlabel("normalized coordinate step  Δθ / range   (g range 1.0, k range 1.8)")
    ax.set_ylabel("α⊥(H4) [s⁻¹]")
    ax.grid(axis="y", color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(loc="upper right", handlelength=1.6, fontsize=5.6, borderaxespad=0.1)
    ax.set_title(
        "D  Boundary motion: port derivative + Newton",
        loc="left",
        color=INK,
    )


def main() -> int:
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.16, 5.8),
        gridspec_kw={"wspace": 0.5, "hspace": 0.42, "width_ratios": [1.0, 1.1]},
    )
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[1, 0])
    panel_d(axes[1, 1])
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        FIG.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={
            "CreationDate": None,
            "ModDate": None,
            "Producer": None,
            "Creator": None,
        },
    )
    fig.savefig(
        FIG.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        metadata={"Software": None},
    )
    fig.savefig(
        FIG.with_suffix(".svg"),
        bbox_inches="tight",
        metadata={"Date": None, "Creator": None},
    )
    print("wrote", FIG.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
