"""Classification metrics against hand-computed cases."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.uncertainty.classification import (
    average_precision,
    confusion,
    permutation_auc,
    precision_at_k,
    recall_at_k,
    roc_auc,
)


def test_perfect_and_reversed_ranking():
    labels = [0, 0, 1, 1]
    assert roc_auc(labels, [0.1, 0.2, 0.8, 0.9]) == pytest.approx(1.0)
    assert roc_auc(labels, [0.9, 0.8, 0.2, 0.1]) == pytest.approx(0.0)


def test_auc_of_a_constant_score_is_one_half():
    assert roc_auc([0, 1, 0, 1], [3.0, 3.0, 3.0, 3.0]) == pytest.approx(0.5)


def test_auc_matches_the_rank_sum_by_hand():
    # scores 1,2,3,4 with positives at 2 and 4: pairs (2>1),(4>1),(4>3) win, (2<3) loses
    assert roc_auc([0, 1, 0, 1], [1.0, 2.0, 3.0, 4.0]) == pytest.approx(3.0 / 4.0)


def test_auc_and_flipped_auc_sum_to_one():
    rng = np.random.default_rng(4)
    labels = rng.integers(0, 2, size=50)
    scores = rng.normal(size=50)
    assert roc_auc(labels, scores) + roc_auc(labels, -scores) == pytest.approx(1.0)


def test_average_precision_by_hand():
    # ranked 1,0,1: precision at the two hits is 1.0 and 2/3
    value = average_precision([1, 0, 1], [0.9, 0.5, 0.1])
    assert value == pytest.approx((1.0 + 2.0 / 3.0) / 2.0)


def test_average_precision_of_a_perfect_ranking_is_one():
    assert average_precision([1, 1, 0, 0], [4, 3, 2, 1]) == pytest.approx(1.0)


def test_precision_and_recall_at_k():
    labels = [1, 0, 1, 0, 1]
    scores = [5, 4, 3, 2, 1]
    assert precision_at_k(labels, scores, 2) == pytest.approx(0.5)
    assert recall_at_k(labels, scores, 2) == pytest.approx(1.0 / 3.0)
    assert recall_at_k(labels, scores, 5) == pytest.approx(1.0)


def test_confusion_and_balanced_accuracy():
    c = confusion([1, 1, 0, 0, 0], [1, 0, 1, 0, 0])
    assert (c.true_positive, c.false_negative, c.false_positive, c.true_negative) == (1, 1, 1, 2)
    assert c.precision == pytest.approx(0.5)
    assert c.recall == pytest.approx(0.5)
    assert c.specificity == pytest.approx(2.0 / 3.0)
    assert c.balanced_accuracy == pytest.approx(0.5 * (0.5 + 2.0 / 3.0))


def test_permutation_null_is_centred_on_one_half():
    rng = np.random.default_rng(11)
    labels = np.array([0, 1] * 25)
    scores = rng.normal(size=50)
    result = permutation_auc(labels, scores, draws=400, rng=rng)
    assert result["null_mean"] == pytest.approx(0.5, abs=0.03)
    assert 0.0 <= result["p_two_sided"] <= 1.0


def test_permutation_p_is_small_for_a_perfect_predictor():
    rng = np.random.default_rng(2)
    labels = np.array([0] * 20 + [1] * 20)
    scores = np.arange(40, dtype=float)
    assert permutation_auc(labels, scores, draws=400, rng=rng)["p_two_sided"] < 0.01
