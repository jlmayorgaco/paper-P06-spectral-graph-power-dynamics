from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import root

from .dae import RealVector, SemiExplicitDAE, residual

_TERMINAL_OK = "SUCCESS"


@dataclass(frozen=True)
class Equilibrium:
    """A solved operating point with its acceptance evidence.

    A sample is never silently accepted: ``status`` is ``SUCCESS`` only when both
    residual norms are below tolerance. The optimizer success flag is NOT part of
    the test, because a derivative-free solver started at the solution reports
    failure to progress while sitting on a machine-precision residual.
    """

    x: RealVector
    z: RealVector
    theta: dict[str, float]
    norm_f: float
    norm_g: float
    status: str
    iterations: int

    @property
    def ok(self) -> bool:
        return self.status == _TERMINAL_OK


def solve_equilibrium(
    model: SemiExplicitDAE,
    theta: dict[str, float],
    x0: RealVector,
    z0: RealVector,
    *,
    tol: float = 1e-10,
    method: str = "hybr",
) -> Equilibrium:
    """Solve ``f=0, g=0`` from a warm start and certify the residual."""

    n_x = int(model.n_x)

    def stacked(w: NDArray[np.float64]) -> NDArray[np.float64]:
        return residual(model, w[:n_x], w[n_x:], theta).stacked

    guess = np.concatenate([np.asarray(x0, float), np.asarray(z0, float)])
    solution = root(stacked, guess, method=method, tol=tol)
    x = np.asarray(solution.x[:n_x], dtype=np.float64)
    z = np.asarray(solution.x[n_x:], dtype=np.float64)
    res = residual(model, x, z, theta)
    if not np.isfinite(res.norm) or res.norm > tol:
        # A derivative-free solver started exactly at the solution reports
        # "not making good progress" and sets success to False, so the optimizer
        # flag cannot be the acceptance test. Retry once with a different method
        # and then judge on the residual, which is what defines an equilibrium.
        retry = root(stacked, guess, method="lm", tol=tol)
        candidate = residual(model, retry.x[:n_x], retry.x[n_x:], theta)
        if np.isfinite(candidate.norm) and candidate.norm < res.norm:
            solution, res = retry, candidate
            x = np.asarray(retry.x[:n_x], dtype=np.float64)
            z = np.asarray(retry.x[n_x:], dtype=np.float64)
    status = _TERMINAL_OK
    if not np.isfinite(res.norm) or res.norm > tol:
        status = "EQ_RESIDUAL_FAIL"
    return Equilibrium(
        x=x,
        z=z,
        theta=dict(theta),
        norm_f=res.norm_f,
        norm_g=res.norm_g,
        status=status,
        iterations=int(solution.get("nfev", -1)),
    )
