"""Numerical method of steps for the exact linear retarded PLL DDE.

The equation is x' = A0*x + B*C.T*x(t-tau). The matrix-exponential step
integrates the undelayed stiff part exactly and linearly interpolates the
known delayed history. The result is V2 linear evidence, never nonlinear DDE.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent


def matrix(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def vector(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", skiprows=1)


def simulate(A: np.ndarray, B: np.ndarray, C: np.ndarray, tau: float,
             root: complex, mode: np.ndarray, steps_per_delay: int,
             horizon: float = 8.0) -> tuple[float, float, int, list[tuple[float,float]]]:
    n = len(mode)
    h = tau / steps_per_delay
    block = np.zeros((n + 20, n + 20))
    block[:n, :n] = A
    block[:n, n:n+10] = B
    block[n:n+10, n+10:] = np.eye(10)
    transition = expm(block * h)
    phi = transition[:n, :n]
    j0 = transition[:n, n:n+10]
    j1 = transition[:n, n+10:] / h
    g0, g1 = j0 - j1, j1
    # Complex history is a legitimate linear combination of two real
    # solutions and allows a clean growth-rate estimate without peak finding.
    history = np.asarray([mode * np.exp(root * ((j - steps_per_delay) * h))
                          for j in range(steps_per_delay + 1)], dtype=np.complex128)
    ring = np.zeros((steps_per_delay + 1, n), dtype=np.complex128)
    ring[:] = history
    kmax = round(horizon / h)
    x = ring[-1].copy()
    initial_norm = np.linalg.norm(x)
    trace=[(0.0,1.0)]
    stride=max(1,round(0.05/h))
    for k in range(kmax):
        q0 = C.T @ ring[k % (steps_per_delay + 1)]
        q1 = C.T @ ring[(k + 1) % (steps_per_delay + 1)]
        x = phi @ x + g0 @ q0 + g1 @ q1
        ring[k % (steps_per_delay + 1)] = x
        if (k+1)%stride==0 or k==kmax-1:
            trace.append(((k+1)*h,float(np.linalg.norm(x)/initial_norm)))
    final_norm = np.linalg.norm(x)
    measured = np.log(final_norm / initial_norm) / (kmax * h)
    return float(measured), float(final_norm / initial_norm), kmax,trace


def main() -> None:
    design_id = sys.argv[1] if len(sys.argv) > 1 else "pending_41ms"
    source = HERE / "linear_dde" / design_id
    A, B, C = (matrix(source / name) for name in ("A0.csv", "B.csv", "C.csv"))
    assert A.shape[0] == A.shape[1] == B.shape[0] == C.shape[0]
    with (source / "roots.csv").open(newline="") as stream:
        roots = list(csv.DictReader(stream))
    result = []
    traces=[]
    for row in roots:
        factor = float(row["factor"])
        tag = str(factor).replace(".", "p")
        mode = vector(source / f"mode_{tag}_real.csv") + 1j * vector(source / f"mode_{tag}_imag.csv")
        root = complex(float(row["root_real"]), float(row["root_imag"]))
        tau = float(row["tau_ms"]) / 1000
        for points in (80, 160):
            measured, ratio, steps, trace = simulate(A, B, C, tau, root, mode, points)
            result.append(dict(design_id=design_id, factor=factor, tau_ms=tau*1000,
                               steps_per_delay=points, time_steps=steps,
                               predicted_growth_s_inv=root.real,
                               measured_growth_s_inv=measured,
                               growth_rate_error_s_inv=measured-root.real,
                               final_initial_amplitude_ratio=ratio,
                               observed_behavior="GROWTH" if measured > 0 else "DECAY",
                               status="LINEAR_EXACT_DDE_METHOD_OF_STEPS_NUMERICAL"))
            print(result[-1], flush=True)
            if points==160:
                traces.extend(dict(design_id=design_id,factor=factor,t_s=t,relative_amplitude=amp,
                                   status="LINEAR_DDE_V2") for t,amp in trace)
    with (HERE / "TABLE_F8_POSITIVE_DELAY_VALIDATION.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    with (HERE / f"TIME_DOMAIN_V2_{design_id}.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(traces[0]))
        writer.writeheader();writer.writerows(traces)


if __name__ == "__main__":
    main()
