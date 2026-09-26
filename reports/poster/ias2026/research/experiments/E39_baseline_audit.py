"""E39 - overnight 11. Independent re-audit of the lower-order baseline result.

The earlier headline was an AUC near 0.14 for the lower-order reconstructions.
An AUC below one half is a statement about the ORDERING of a predictor on THIS
benchmark, not a proof that a method is universally anti-predictive, and the
wording rule in this experiment enforces that distinction.

Every predictor is scored on its RAW value, and both the AUC and its complement
are reported, so no hidden sign convention can flatter or damn anything. The
conventional expectation for each predictor is stated alongside.

Portfolio sizes 4, 5 and 6 are audited separately, because the base rate moves
from 8 % to 93 % across them and a single pooled number would be meaningless.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _overnight import Experiment, pin_blas_threads
from ibr_cycles.uncertainty.classification import (
    average_precision,
    confusion,
    permutation_auc,
    precision_at_k,
    recall_at_k,
    roc_auc,
)

NAME = "E39_baseline_audit"
SIZES = (4, 5, 6)
PERMUTATIONS = 5000
#: The direction a practitioner would expect the raw value to move with risk.
EXPECTATION = {
    "additive_prediction": "higher is more unstable",
    "pairwise_prediction": "higher is more unstable",
    "replaced_mw": "higher is more unstable",
    "replaced_inertia_fraction": "higher is more unstable",
    "min_scr": "LOWER is more unstable",
    "mean_scr": "LOWER is more unstable",
    "gscr": "LOWER is more unstable",
    "max_miif": "higher is more unstable",
}


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Do the conventional and lower-order baselines rank instability?",
        config={
            "predictors": list(EXPECTATION),
            "expectation": EXPECTATION,
            "sizes": list(SIZES),
            "permutations": PERMUTATIONS,
            "threshold_rule": "lower-order reconstruction predicts unstable iff alpha > 0",
            "wording_rule": (
                "an AUC below one half with a complement above one half is reported as "
                "'inversely ordered on this benchmark', never as 'universally "
                "anti-predictive'"
            ),
        },
        workers=1,
    )
    started = time.time()
    rng = np.random.default_rng(20260917)
    data = pd.read_csv(RESULTS / "tables" / "E13_baseline_challenge_predictors.csv")

    rows = []
    for size in SIZES:
        block = data[data["size"] == size]
        labels = block.unstable.to_numpy(bool)
        if labels.sum() in (0, labels.size):
            experiment.note(f"size {size}: degenerate labels, skipped")
            continue
        for predictor in EXPECTATION:
            scores = block[predictor].to_numpy(float)
            auc = roc_auc(labels, scores)
            permutation = permutation_auc(labels, scores, draws=PERMUTATIONS, rng=rng)
            row = {
                "size": size,
                "n": int(labels.size),
                "unstable": int(labels.sum()),
                "base_rate": float(labels.mean()),
                "predictor": predictor,
                "expectation": EXPECTATION[predictor],
                "roc_auc": auc,
                "roc_auc_sign_flipped": 1.0 - auc,
                "pr_auc": average_precision(labels, scores),
                "pr_auc_sign_flipped": average_precision(labels, -scores),
                "precision_at_10": precision_at_k(labels, scores, 10),
                "recall_at_10": recall_at_k(labels, scores, 10),
                "precision_at_10_flipped": precision_at_k(labels, -scores, 10),
                "recall_at_10_flipped": recall_at_k(labels, -scores, 10),
                "permutation_null_mean": permutation["null_mean"],
                "permutation_p": permutation["p_two_sided"],
            }
            if predictor in ("additive_prediction", "pairwise_prediction"):
                predicted = scores > 0.0
                matrix = confusion(labels, predicted)
                row.update(
                    {f"threshold_{k}": v for k, v in matrix.as_dict().items()}
                )
                row["threshold_rule"] = "predicted unstable iff reconstructed alpha > 0"
            rows.append(row)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E39_baseline_audit.csv")

    lower_order = table[
        table.predictor.isin(("additive_prediction", "pairwise_prediction"))
    ]
    conventional = table[
        table.predictor.isin(("min_scr", "mean_scr", "gscr", "max_miif"))
    ]
    inverted = table[(table.roc_auc < 0.5) & (table.permutation_p < 0.05)]
    checks = {
        "rows": int(len(table)),
        "lower_order_auc_range": [
            float(lower_order.roc_auc.min()),
            float(lower_order.roc_auc.max()),
        ],
        "lower_order_flipped_auc_range": [
            float(lower_order.roc_auc_sign_flipped.min()),
            float(lower_order.roc_auc_sign_flipped.max()),
        ],
        "conventional_auc_range": [
            float(conventional.roc_auc.min()),
            float(conventional.roc_auc.max()),
        ],
        "predictors_significantly_inverted": int(len(inverted)),
        "inverted_list": sorted(
            {f"{r.predictor}@{r['size']}" for _, r in inverted.iterrows()}
        ),
        "permutation_null_mean_overall": float(table.permutation_null_mean.mean()),
        "any_predictor_above_0_8": sorted(
            {
                f"{r.predictor}@{r['size']}"
                for _, r in table.iterrows()
                if max(r.roc_auc, r.roc_auc_sign_flipped) > 0.8
            }
        ),
    }

    print()
    print(
        table[
            [
                "size",
                "predictor",
                "base_rate",
                "roc_auc",
                "roc_auc_sign_flipped",
                "pr_auc",
                "precision_at_10",
                "permutation_p",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:8.4f}")
    )
    print()
    thresholds = table[table.threshold_rule.notna()] if "threshold_rule" in table else pd.DataFrame()
    if len(thresholds):
        print(
            thresholds[
                [
                    "size",
                    "predictor",
                    "threshold_tp",
                    "threshold_fp",
                    "threshold_tn",
                    "threshold_fn",
                    "threshold_precision",
                    "threshold_recall",
                    "threshold_balanced_accuracy",
                ]
            ].to_string(index=False, float_format=lambda v: f"{v:8.4f}")
        )
    print()
    for key, value in checks.items():
        print(f"  {key:36s} {value}")
    experiment.finish("REPORTED", checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
