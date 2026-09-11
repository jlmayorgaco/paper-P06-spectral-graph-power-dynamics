"""UC02 figure: E39 ROC/PR re-drawn with the unit-correct predictors (new version).

The frozen E39_ROC_PR.png labels the converter rating (MVA) as "replaced_mw".
This figure shows it as "retired rating Sn [MVA]" next to the measured
"removed dispatch Pg [MW]". Source data are written beside the figure.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from _uc import FROZEN_TABLES, RESULTS, out_dir  # noqa: E402

OUT = out_dir("UC02")
SERIES = (  # fixed categorical order; identity by label, never by rank
    ("replaced_inertia_fraction", "removed inertia fraction", "#1b9e77", "-"),
    ("replaced_sn_mva", "retired rating Sn [MVA] (old 'replaced MW')", "#1f5fa8", "-"),
    ("replaced_pg_mw", "removed dispatch Pg [MW]", "#d95f02", "-"),
    ("min_scr", "min SCR (raw)", "#7a3fb0", "--"),
    ("gscr", "gSCR (raw)", "#6b6b6b", "--"),
    ("additive_prediction", "additive reconstruction", "#c2185b", ":"),
    ("pairwise_prediction", "pairwise reconstruction", "#e6a100", ":"),
)


def roc(labels, scores):
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    tpr = np.concatenate([[0], np.cumsum(y) / y.sum()])
    fpr = np.concatenate([[0], np.cumsum(~y) / (~y).sum()])
    return fpr, tpr


def pr(labels, scores):
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    tp = np.cumsum(y)
    return tp / y.sum(), tp / np.arange(1, y.size + 1)


def main() -> int:
    data = pd.read_csv(FROZEN_TABLES / "E13_baseline_challenge_predictors.csv")
    q = pd.read_csv(RESULTS / "UC01" / "UC01_census_quantities.csv")
    data = data.merge(q[["members", "replaced_sn_mva", "replaced_pg_mw"]], on="members")
    scores = pd.read_csv(OUT / "UC02_baseline_scores.csv")
    auc = {(r["size"], r.predictor): r.roc_auc for _, r in scores.iterrows()}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6), constrained_layout=True)
    source = []
    for j, size in enumerate((4, 5, 6)):
        block = data[data["size"] == size]
        labels = block.unstable.to_numpy(bool)
        for col, label, color, style in SERIES:
            s = block[col].to_numpy(float)
            fpr, tpr = roc(labels, s)
            rec, prec = pr(labels, s)
            key = "replaced_mw" if col == "replaced_sn_mva" else col
            value = auc.get((size, col), auc.get((size, key)))
            axes[0, j].plot(
                fpr, tpr, color=color, ls=style, lw=2, label=f"{label} ({value:.2f})"
            )
            axes[1, j].plot(rec, prec, color=color, ls=style, lw=2)
            source += [
                {"size": size, "predictor": col, "panel": "roc", "x": a, "y": b}
                for a, b in zip(fpr, tpr, strict=True)
            ]
            source += [
                {"size": size, "predictor": col, "panel": "pr", "x": a, "y": b}
                for a, b in zip(rec, prec, strict=True)
            ]
        rate = labels.mean()
        axes[0, j].plot([0, 1], [0, 1], color="#999999", lw=1, ls="--")
        axes[1, j].axhline(rate, color="#999999", lw=1, ls="--")
        axes[0, j].set_title(f"ROC, size {size} (base rate {rate:.2f})")
        axes[1, j].set_title(f"precision-recall, size {size}")
        axes[0, j].set_xlabel("false positive rate")
        axes[0, j].set_ylabel("true positive rate")
        axes[1, j].set_xlabel("recall")
        axes[1, j].set_ylabel("precision")
        axes[0, j].legend(fontsize=7, loc="lower right")
        for ax in axes[:, j]:
            ax.grid(alpha=0.3)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
    fig.suptitle(
        "Baseline audit, unit-corrected (UC02): rating Sn [MVA] and dispatch Pg [MW] "
        "are different predictors; AUC in parentheses (raw value)"
    )
    fig.savefig(OUT / "UC02_ROC_PR_unit_corrected.png", dpi=150)
    pd.DataFrame(source).to_csv(OUT / "UC02_figure_source.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
