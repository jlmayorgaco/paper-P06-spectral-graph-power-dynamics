"""NL02: derivatives of the algebraic manifold and of the reduced field (v2 N.47-N.50).

In a fixed active set with M = I (L0), let z = psi(x) solve g(x, z) = 0 and
r(x) = f(x, psi(x)). With G_z factorized once (LU), never inverted:

    J_psi = -G_z^-1 G_x                                  first order of the chart
    A     = f_x + f_z J_psi                              reduced Jacobian
    H_psi[a, b] = -G_z^-1 D2g[(a, J_psi a), (b, J_psi b)]
    H_r[a, b]   = D2f[(a, J_psi a), (b, J_psi b)] + f_z H_psi[a, b]

Inputs u and outputs h(x, z, u):

    B = f_u + f_z (-G_z^-1 G_u),  C = h_x + h_z J_psi,  D = h_u + h_z (-G_z^-1 G_u).

Second derivatives are directional (Hessian-vector products) from a fourth-point
central difference of the REAL residuals:

    D2F[p, q] ~ [F(y + h p + h q) - F(y + h p - h q) - F(y - h p + h q)
                 + F(y - h p - h q)] / (4 h^2),        y = (x, z).

No complex step: the model contains conjugates and moduli. No dense cubic
tensor is formed. Verification compares H_r with a direct second difference of
r, in which psi(x) is recomputed by Newton at every evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from ..dynamics.linearize import central_difference_jacobians


@dataclass
class Chart:
    x0: np.ndarray
    z0: np.ndarray
    fx: np.ndarray
    fz: np.ndarray
    gx: np.ndarray
    gz: np.ndarray
    lu: tuple

    @property
    def j_psi(self) -> np.ndarray:
        return -lu_solve(self.lu, self.gx)

    @property
    def a(self) -> np.ndarray:
        return self.fx + self.fz @ self.j_psi


def chart(dae, x0, z0) -> Chart:
    jac = central_difference_jacobians(dae, x0, z0, {})
    return Chart(
        np.asarray(x0, float),
        np.asarray(z0, float),
        jac.fx,
        jac.fz,
        jac.gx,
        jac.gz,
        lu_factor(jac.gz),
    )


def second_directional(fun, y, p, q, h):
    return (
        fun(y + h * p + h * q)
        - fun(y + h * p - h * q)
        - fun(y - h * p + h * q)
        + fun(y - h * p - h * q)
    ) / (4.0 * h * h)


def lifted(ch: Chart, a: np.ndarray) -> np.ndarray:
    """(a, J_psi a): a state direction lifted to the tangent of the manifold."""

    return np.concatenate([a, ch.j_psi @ a])


def hessian_reduced(
    dae, ch: Chart, a, b, h: float = 1e-4
) -> tuple[np.ndarray, np.ndarray]:
    """(H_psi[a, b], H_r[a, b]) by the manifold formulas."""

    n_x = ch.x0.size
    y0 = np.concatenate([ch.x0, ch.z0])
    pa, pb = lifted(ch, a), lifted(ch, b)

    def f_of(y):
        return dae.f(y[:n_x], y[n_x:], {})

    def g_of(y):
        return dae.g(y[:n_x], y[n_x:], {})

    d2g = second_directional(g_of, y0, pa, pb, h)
    d2f = second_directional(f_of, y0, pa, pb, h)
    h_psi = -lu_solve(ch.lu, d2g)
    h_r = d2f + ch.fz @ h_psi
    return h_psi, h_r


def psi(dae, x, z_start, lu, tol=1e-12, max_iter=60) -> np.ndarray:
    """Solve g(x, z) = 0 near z_start: chord iterations on the frozen LU, then Newton.

    Converged when the residual is below ``tol`` or the update has stalled at
    rounding level (|dz| <= 1e-15 (1 + |z|)) with residual below 1e-9. Otherwise
    the point is declared outside the chart (RuntimeError); it is never
    returned silently.
    """

    z = np.array(z_start, dtype=float)

    def done(res, dz):
        r = np.abs(res).max()
        return r < tol or (
            np.abs(dz).max() <= 1e-15 * (1 + np.abs(z).max()) and r < 1e-9
        )

    for _ in range(max_iter):
        res = dae.g(x, z, {})
        dz = lu_solve(lu, res)
        if done(res, dz):
            return z
        z = z - dz
    jac = central_difference_jacobians(dae, x, z, {})
    for _ in range(20):
        res = dae.g(x, z, {})
        dz = np.linalg.solve(jac.gz, res)
        if done(res, dz):
            return z
        z = z - dz
    raise RuntimeError(
        "algebraic chart: no convergence (singular or outside the chart)"
    )


def reduced_field(dae, ch: Chart, x) -> np.ndarray:
    z = psi(dae, x, ch.z0, ch.lu)
    return dae.f(x, z, {})


def hessian_direct(dae, ch: Chart, a, b, h: float) -> np.ndarray:
    """Second directional difference of r(x) = f(x, psi(x)), psi solved each time."""

    return second_directional(lambda x: reduced_field(dae, ch, x), ch.x0, a, b, h)


def input_output(model, ch: Chart, du: float = 1e-6):
    """(B, C, D) for the declared inputs and outputs of a PhasorModel."""

    x0, z0 = ch.x0, ch.z0
    n_u = model.n_u
    f_u = np.zeros((x0.size, n_u))
    g_u = np.zeros((z0.size, n_u))
    for j in range(n_u):
        e = np.zeros(n_u)
        e[j] = du
        f_u[:, j] = (
            model.residual_f(x0, z0, e).values - model.residual_f(x0, z0, -e).values
        ) / (2 * du)
        g_u[:, j] = (
            model.residual_g(x0, z0, e).values - model.residual_g(x0, z0, -e).values
        ) / (2 * du)
    h0 = model.outputs(x0, z0)
    n_h = h0.values.size
    h_x = np.zeros((n_h, x0.size))
    h_z = np.zeros((n_h, z0.size))
    h_u = np.zeros((n_h, n_u))
    step = 1e-6
    for j in range(x0.size):
        e = np.zeros(x0.size)
        e[j] = step
        h_x[:, j] = (
            model.outputs(x0 + e, z0).values - model.outputs(x0 - e, z0).values
        ) / (2 * step)
    for j in range(z0.size):
        e = np.zeros(z0.size)
        e[j] = step
        h_z[:, j] = (
            model.outputs(x0, z0 + e).values - model.outputs(x0, z0 - e).values
        ) / (2 * step)
    for j in range(n_u):
        e = np.zeros(n_u)
        e[j] = du
        h_u[:, j] = (
            model.outputs(x0, z0, e).values - model.outputs(x0, z0, -e).values
        ) / (2 * du)
    psi_u = -lu_solve(ch.lu, g_u)
    b = f_u + ch.fz @ psi_u
    c = h_x + h_z @ ch.j_psi
    d = h_u + h_z @ psi_u
    return {"B": b, "C": c, "D": d, "names_u": model.input_names, "names_h": h0.names}
