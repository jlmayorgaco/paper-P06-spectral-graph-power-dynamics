"""Ranking and classification metrics, implemented rather than imported.

scikit-learn is not in the pinned environment and nothing may be installed, so
the handful of metrics the baseline audit needs are written here and checked in
``tests/test_classification.py`` against hand-computed cases.

A note on direction. A predictor that ranks the wrong way round produces an AUC
below one half. That is a statement about THIS benchmark's ordering, not a proof
that the predictor is universally anti-predictive, so ``roc_auc`` is always
reported beside ``1 - roc_auc`` and the caller is expected to say which it means.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


def _ranks(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Average ranks, so ties do not bias the statistic."""

    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    sorted_values = values[order]
    i = 0
    while i < values.size:
        j = i
        while j + 1 < values.size and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def roc_auc(labels, scores) -> float:
    """Area under the ROC curve, by the rank-sum identity, ties averaged."""

    y = np.asarray(labels, dtype=bool)
    s = np.asarray(scores, dtype=float)
    positives = int(y.sum())
    negatives = int((~y).sum())
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = _ranks(s)
    return float((ranks[y].sum() - positives * (positives + 1) / 2.0) / (positives * negatives))


def average_precision(labels, scores) -> float:
    """Area under the precision-recall curve, the step-wise average precision."""

    y = np.asarray(labels, dtype=bool)
    s = np.asarray(scores, dtype=float)
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    cumulative = np.cumsum(y)
    precision = cumulative / np.arange(1, y.size + 1)
    return float(precision[y].sum() / y.sum())


def precision_at_k(labels, scores, k: int) -> float:
    y = np.asarray(labels, dtype=bool)
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    k = min(k, y.size)
    return float(y[order][:k].sum() / k) if k else float("nan")


def recall_at_k(labels, scores, k: int) -> float:
    y = np.asarray(labels, dtype=bool)
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    return float(y[order][: min(k, y.size)].sum() / y.sum())


@dataclass(frozen=True)
class Confusion:
    true_positive: int
    false_positive: int
    true_negative: int
    false_negative: int

    @property
    def precision(self) -> float:
        denominator = self.true_positive + self.false_positive
        return self.true_positive / denominator if denominator else float("nan")

    @property
    def recall(self) -> float:
        denominator = self.true_positive + self.false_negative
        return self.true_positive / denominator if denominator else float("nan")

    @property
    def specificity(self) -> float:
        denominator = self.true_negative + self.false_positive
        return self.true_negative / denominator if denominator else float("nan")

    @property
    def balanced_accuracy(self) -> float:
        return 0.5 * (self.recall + self.specificity)

    def as_dict(self) -> dict[str, float]:
        return {
            "tp": self.true_positive,
            "fp": self.false_positive,
            "tn": self.true_negative,
            "fn": self.false_negative,
            "precision": self.precision,
            "recall": self.recall,
            "specificity": self.specificity,
            "balanced_accuracy": self.balanced_accuracy,
        }


def confusion(labels, predictions) -> Confusion:
    y = np.asarray(labels, dtype=bool)
    p = np.asarray(predictions, dtype=bool)
    return Confusion(
        true_positive=int((y & p).sum()),
        false_positive=int((~y & p).sum()),
        true_negative=int((~y & ~p).sum()),
        false_negative=int((y & ~p).sum()),
    )


def permutation_auc(labels, scores, *, draws: int, rng) -> dict[str, float]:
    """Null distribution of the AUC under label permutation."""

    y = np.asarray(labels, dtype=bool)
    observed = roc_auc(y, scores)
    null = np.empty(draws)
    for i in range(draws):
        null[i] = roc_auc(rng.permutation(y), scores)
    finite = null[np.isfinite(null)]
    return {
        "observed": observed,
        "null_mean": float(finite.mean()),
        "null_std": float(finite.std()),
        "p_two_sided": float(
            (np.abs(finite - 0.5) >= abs(observed - 0.5)).mean()
        ),
    }
