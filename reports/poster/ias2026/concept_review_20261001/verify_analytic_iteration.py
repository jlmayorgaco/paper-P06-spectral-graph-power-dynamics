"""Check derivation identities, not power-system performance or robustness.

Checks nonlinear index-one DAE sensitivities against central differences of
separately integrated trajectories, including parameter-dependent equilibrium.
Also checks the closed active-set step against the full saddle-point solve.
The small nonlinear DAE is an algebraic verification example, not IEEE-39.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def algebraic(theta, rho, w):
    # Unique exact real solution of z + 0.2*z**3 = b.
    b = 0.15 * theta + 0.2 * rho + w
    return 2 * np.sqrt(5 / 3) * np.sinh(np.arcsinh(1.5 * b * np.sqrt(3 / 5)) / 3)


def terms(t, x, p):
    rho, kp, ki = p
    angle, omega, xi = x
    w = 0.04 * np.sin(1.1 * t) + 0.02 * (1 - np.exp(-t))
    z = algebraic(angle, rho, w)
    mass = 1 - 0.6 * rho
    error = z - np.sin(angle)
    numerator = xi + kp * error - omega
    f = np.array([omega, numerator / mass, ki * error])
    fx = np.array([[0, 1, 0], [-kp * np.cos(angle) / mass, -1 / mass, 1 / mass],
                   [-ki * np.cos(angle), 0, 0]])
    fz = np.array([0, kp / mass, ki])
    fp = np.array([[0, 0, 0], [0.6 * numerator / mass**2, error / mass, 0], [0, 0, error]])
    gx = np.array([-0.15, 0, 0])
    gp = np.array([-0.2, 0, 0])
    gz = 1 + 0.6 * z**2
    return f, fx, fz, fp, gx, gp, gz, z


def initial(p):
    angle = brentq(lambda a: np.sin(a) + 0.2 * np.sin(a)**3 - 0.15 * a - 0.2 * p[0],
                   -0.5, 0.5, xtol=1e-14)
    x = np.array([angle, 0.0, 0.0])
    f, fx, fz, fp, gx, gp, gz, z = terms(0.0, x, p)
    eq_jac = np.block([[fx, fz[:, None]], [gx[None, :], np.array([[gz]])]])
    sens = np.linalg.solve(eq_jac, -np.vstack([fp, gp]))
    return x, sens[:3], float(np.max(np.abs(f)))


def run(p, sensitivities):
    x0, s0, residual = initial(p)
    def rhs(t, y):
        f, fx, fz, fp, gx, gp, gz, z = terms(t, y[:3], p)
        if not sensitivities:
            return f
        s = y[3:].reshape(3, 3)
        ds = (fx - np.outer(fz, gx) / gz) @ s + fp - np.outer(fz, gp) / gz
        return np.r_[f, ds.ravel()]
    y0 = np.r_[x0, s0.ravel()] if sensitivities else x0
    sol = solve_ivp(rhs, (0, 4.0), y0, method="DOP853", rtol=2e-12, atol=2e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    x = sol.y[:3, -1]
    *_, gx, gp, gz, z = terms(4.0, x, p)
    out = np.r_[x, z]
    if not sensitivities:
        return out
    s = sol.y[3:, -1].reshape(3, 3)
    zs = -(gx @ s + gp) / gz
    return out, np.vstack([s, zs]), residual


def main():
    p = np.array([0.4, 2.0, 0.7])
    out, derivative, initial_residual = run(p, True)
    comparisons = []
    for eps in (2e-4, 1e-4, 5e-5):
        fd = np.column_stack([(run(p + eps * e, False) - run(p - eps * e, False)) / (2 * eps)
                              for e in np.eye(3)])
        comparisons.append({"step": eps, "max_absolute_error": float(np.max(np.abs(fd - derivative))),
                            "relative_frobenius_error": float(np.linalg.norm(fd - derivative) / np.linalg.norm(derivative))})

    rng = np.random.default_rng(701)
    n, m = 6, 3
    r = rng.normal(size=(n, n))
    h = r.T @ r + np.eye(n)
    a = rng.normal(size=(m, n))
    c = rng.normal(size=n)
    multiplier = 0.5 + rng.random(m)
    target_step = np.linalg.solve(h, c - a.T @ multiplier)
    g = -a @ target_step
    hinv_c = np.linalg.solve(h, c)
    hinv_at = np.linalg.solve(h, a.T)
    mu = np.linalg.solve(a @ hinv_at, g + a @ hinv_c)
    step = hinv_c - hinv_at @ mu
    saddle = np.block([[h, a.T], [a, np.zeros((m, m))]])
    reference = np.linalg.solve(saddle, np.r_[c, -g])
    qp = {"step_error_vs_saddle_solve": float(np.max(np.abs(step - reference[:n]))),
          "stationarity_residual": float(np.max(np.abs(h @ step - c + a.T @ mu))),
          "active_constraint_residual": float(np.max(np.abs(g + a @ step))),
          "min_multiplier": float(mu.min())}
    result = {"scope": "DERIVATION_IDENTITY_CHECK_ONLY; not a physical SG-GFL experiment or robust certificate",
              "initial_equilibrium_residual": initial_residual,
              "nonlinear_dae_gradient_checks": comparisons,
              "active_set_formula_check": qp}
    assert comparisons[-1]["relative_frobenius_error"] < 2e-6
    assert qp["step_error_vs_saddle_solve"] < 1e-11
    assert qp["stationarity_residual"] < 1e-11
    assert qp["active_constraint_residual"] < 1e-11
    assert qp["min_multiplier"] > 0
    destination = Path(__file__).with_name("ANALYTIC_ITERATION_CHECK.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
