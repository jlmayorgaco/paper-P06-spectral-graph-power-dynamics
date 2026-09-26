"""Regression estimators checked against closed-form and known answers."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import t as student_t

from ibr_cycles.uncertainty.regression import least_squares, logistic

NAMES = ("const", "x")


def simple(x, y):
    return least_squares(np.column_stack([np.ones_like(x), x]), y, NAMES)


def test_least_squares_matches_the_closed_form_simple_regression():
    rng = np.random.default_rng(7)
    x = rng.normal(size=60)
    y = 1.5 - 2.25 * x + 0.3 * rng.normal(size=60)
    fit = simple(x, y)

    slope = float(np.cov(x, y, ddof=1)[0, 1] / np.var(x, ddof=1))
    intercept = float(y.mean() - slope * x.mean())
    assert fit.coefficients[1] == pytest.approx(slope, rel=1e-12)
    assert fit.coefficients[0] == pytest.approx(intercept, rel=1e-12)

    residual = y - intercept - slope * x
    variance = float(residual @ residual) / (x.size - 2)
    error = np.sqrt(variance / float(np.sum((x - x.mean()) ** 2)))
    assert fit.standard_errors[1] == pytest.approx(error, rel=1e-12)
    assert fit.p_values[1] == pytest.approx(
        2.0 * student_t.sf(abs(slope / error), x.size - 2), rel=1e-12
    )


def test_least_squares_is_exact_on_noise_free_data():
    x = np.linspace(-1.0, 1.0, 25)
    fit = simple(x, 4.0 - 3.0 * x)
    assert fit.coefficients == pytest.approx([4.0, -3.0], abs=1e-12)
    assert fit.goodness == pytest.approx(1.0, abs=1e-12)


def test_least_squares_p_value_is_large_for_a_pure_null():
    rng = np.random.default_rng(11)
    x = rng.normal(size=200)
    fit = simple(x, rng.normal(size=200))
    assert fit.p_of("x") > 0.05


def test_logistic_recovers_the_generating_coefficients():
    rng = np.random.default_rng(3)
    x = rng.normal(size=4000)
    truth = (0.4, 1.8)
    probability = 1.0 / (1.0 + np.exp(-(truth[0] + truth[1] * x)))
    y = (rng.uniform(size=x.size) < probability).astype(float)
    fit = logistic(np.column_stack([np.ones_like(x), x]), y, NAMES)

    assert fit.converged
    assert fit.coefficients[0] == pytest.approx(truth[0], abs=0.12)
    assert fit.coefficients[1] == pytest.approx(truth[1], abs=0.20)
    assert fit.p_of("x") < 1e-20


def test_logistic_score_vanishes_at_the_reported_solution():
    rng = np.random.default_rng(5)
    x = rng.normal(size=300)
    y = (rng.uniform(size=300) < 1.0 / (1.0 + np.exp(-x))).astype(float)
    design = np.column_stack([np.ones_like(x), x])
    fit = logistic(design, y, NAMES)

    probability = 1.0 / (1.0 + np.exp(-(design @ fit.coefficients)))
    assert np.max(np.abs(design.T @ (y - probability))) < 1e-8


def test_logistic_standard_error_matches_the_observed_information():
    rng = np.random.default_rng(13)
    x = rng.normal(size=500)
    y = (rng.uniform(size=500) < 1.0 / (1.0 + np.exp(-(0.5 * x)))).astype(float)
    design = np.column_stack([np.ones_like(x), x])
    fit = logistic(design, y, NAMES)

    probability = 1.0 / (1.0 + np.exp(-(design @ fit.coefficients)))
    weight = probability * (1.0 - probability)
    information = design.T @ (design * weight[:, None])
    expected = np.sqrt(np.diag(np.linalg.inv(information)))
    assert fit.standard_errors == pytest.approx(expected, rel=1e-5)


def test_logistic_survives_a_separable_design_without_blowing_up():
    x = np.linspace(-1.0, 1.0, 40)
    y = (x > 0).astype(float)
    fit = logistic(np.column_stack([np.ones_like(x), x]), y, NAMES)
    assert np.all(np.isfinite(fit.coefficients))
    assert fit.coefficients[1] > 0.0
