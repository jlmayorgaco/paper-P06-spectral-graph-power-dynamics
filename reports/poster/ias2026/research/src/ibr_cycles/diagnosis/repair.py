from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize


@dataclass(frozen=True)
class RepairResult:
    """Outcome of a minimum-change repair of a destabilized portfolio."""

    strategy: str
    theta: dict[str, float]
    change_norm: float
    spectral_abscissa: float
    damping: float
    assets_retained: tuple[str, ...]
    success: bool
    message: str

    @property
    def meets_target(self) -> bool:
        return self.success and self.spectral_abscissa < 0.0


def scalar_repair(
    abscissa: Callable[[float], float],
    bracket: tuple[float, float],
    *,
    target: float = -0.02,
    name: str = "theta",
    tolerance: float = 1e-12,
) -> tuple[float, float]:
    """Find the critical value and the value achieving a target abscissa.

    Returns (critical, repaired). Both are bisected on a bracket the caller must
    verify brackets a sign change; no extrapolation is performed.
    """

    critical = float(brentq(abscissa, *bracket, xtol=tolerance))
    repaired = float(brentq(lambda v: abscissa(v) - target, *bracket, xtol=tolerance))
    return critical, repaired


def minimum_norm_repair(
    abscissa: Callable[[NDArray[np.float64]], float],
    names: Sequence[str],
    baseline: NDArray[np.float64],
    *,
    target: float = -0.02,
    weights: NDArray[np.float64] | None = None,
    bounds: Sequence[tuple[float, float]] | None = None,
    assets_retained: Sequence[str] = (),
    strategy: str = "core_retune",
) -> RepairResult:
    """Minimize the weighted parameter change subject to a stability margin.

    The comparison of interest is not whether a repair exists but whether a
    repair confined to the identified feedback core beats removing a useful
    asset or retuning globally, at equal or smaller change norm.
    """

    w = np.ones_like(baseline) if weights is None else np.asarray(weights, float)

    def objective(theta: NDArray[np.float64]) -> float:
        delta = (theta - baseline) * w
        return float(delta @ delta)

    constraints = [{"type": "ineq", "fun": lambda theta: target - abscissa(theta)}]
    solution = minimize(
        objective,
        baseline.copy(),
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 400, "ftol": 1e-12},
    )
    theta = np.asarray(solution.x, dtype=np.float64)
    achieved = abscissa(theta)
    return RepairResult(
        strategy=strategy,
        theta={name: float(value) for name, value in zip(names, theta, strict=True)},
        change_norm=float(np.linalg.norm((theta - baseline) * w)),
        spectral_abscissa=float(achieved),
        damping=float("nan"),
        assets_retained=tuple(assets_retained),
        success=bool(solution.success and achieved <= target + 1e-9),
        message=str(solution.message),
    )
