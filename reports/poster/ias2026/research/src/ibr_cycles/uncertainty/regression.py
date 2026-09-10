"""Ordinary least squares and logistic regression with standard errors.

Two textbook estimators, implemented here rather than pulled in as a new
dependency in the middle of a preregistered campaign. Both are covered by tests
against closed-form cases in ``tests/test_regression.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.stats import norm, t as student_t


@dataclass(frozen=True)
class Fit:
    """Coefficients with standard errors and two-sided p-values."""

    names: tuple[str, ...]
    coefficients: NDArray[np.float64]
    standard_errors: NDArray[np.float64]
    p_values: NDArray[np.float64]
    n: int
    goodness: float
    converged: bool = True

    def summary(self) -> dict[str, dict[str, float]]:
        return {
            name: {
                "coefficient": float(c),
                "standard_error": float(s),
                "p_value": float(p),
            }
            for name, c, s, p in zip(
                self.names,
                self.coefficients,
                self.standard_errors,
                self.p_values,
                strict=True,
            )
        }

    def index(self, name: str) -> int:
        return self.names.index(name)

    def p_of(self, name: str) -> float:
        return float(self.p_values[self.index(name)])


def least_squares(
    design: NDArray[np.float64], response: NDArray[np.float64], names: tuple[str, ...]
) -> Fit:
    """Ordinary least squares with t-based p-values."""

    x = np.asarray(design, dtype=np.float64)
    y = np.asarray(response, dtype=np.float64)
    n, p = x.shape
    gram_inverse = np.linalg.pinv(x.T @ x)
    beta = gram_inverse @ x.T @ y
    residual = y - x @ beta
    degrees = max(n - p, 1)
    variance = float(residual @ residual) / degrees
    errors = np.sqrt(np.maximum(np.diag(variance * gram_inverse), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = np.where(errors > 0, beta / errors, 0.0)
    p_values = 2.0 * student_t.sf(np.abs(statistic), degrees)
    centred = y - y.mean()
    total = float(centred @ centred)
    r_squared = 1.0 - float(residual @ residual) / total if total > 0 else float("nan")
    return Fit(
        names=names,
        coefficients=beta,
        standard_errors=errors,
        p_values=p_values,
        n=n,
        goodness=r_squared,
    )


def logistic(
    design: NDArray[np.float64],
    response: NDArray[np.float64],
    names: tuple[str, ...],
    *,
    iterations: int = 60,
    ridge: float = 1e-8,
) -> Fit:
    """Logistic regression by iteratively reweighted least squares.

    A small ridge term keeps the Fisher information invertible when the design is
    nearly collinear, which a quadratic surface in two correlated variables often
    is. Wald p-values come from the inverse information at the solution.
    """

    x = np.asarray(design, dtype=np.float64)
    y = np.asarray(response, dtype=np.float64)
    n, p = x.shape
    beta = np.zeros(p)
    converged = False
    information = np.eye(p)
    for _ in range(iterations):
        eta = np.clip(x @ beta, -30.0, 30.0)
        probability = 1.0 / (1.0 + np.exp(-eta))
        weight = np.clip(probability * (1.0 - probability), 1e-9, None)
        information = x.T @ (x * weight[:, None]) + ridge * np.eye(p)
        gradient = x.T @ (y - probability)
        step = np.linalg.solve(information, gradient)
        beta = beta + step
        if float(np.max(np.abs(step))) < 1e-10:
            converged = True
            break
    covariance = np.linalg.pinv(information)
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = np.where(errors > 0, beta / errors, 0.0)
    p_values = 2.0 * norm.sf(np.abs(statistic))
    eta = np.clip(x @ beta, -30.0, 30.0)
    probability = 1.0 / (1.0 + np.exp(-eta))
    epsilon = 1e-12
    log_likelihood = float(
        y @ np.log(probability + epsilon) + (1 - y) @ np.log(1 - probability + epsilon)
    )
    rate = float(np.mean(y))
    null = (
        float(n * (rate * np.log(rate + epsilon) + (1 - rate) * np.log(1 - rate + epsilon)))
        if 0.0 < rate < 1.0
        else float("nan")
    )
    pseudo = 1.0 - log_likelihood / null if np.isfinite(null) and null != 0 else float("nan")
    return Fit(
        names=names,
        coefficients=beta,
        standard_errors=errors,
        p_values=p_values,
        n=n,
        goodness=pseudo,
        converged=converged,
    )
