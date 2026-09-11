"""FC07 (campaign step 10): second-order DAE curvature audit, linear vs full vs no-KCL.

Same frozen DAE. On the algebraic manifold z = psi(x, u):

    D psi       = -G_z^-1 [G_x, G_u]
    D^2 psi[a,b] = -G_z^-1 D^2G[Xi a, Xi b],   Xi a = (a_x, D psi a, a_u)
    D^2 r[a,b]   = D^2F[Xi a, Xi b] + F_z D^2 psi[a,b]

Three reduced models of the response to the same input u(t) (D1 load pulse at
bus 16, 0.5 s, amplitude eps):

    M0  x1' = A x1 + B u                                   (linear)
    M1  x2' = A x2 + 1/2 D^2 r[zeta1, zeta1],  zeta1 = (x1, u)   (full second order)
    M2  the same with F_z D^2 psi removed (algebraic / KCL curvature dropped)

The outputs are evaluated on the reconstructed (x, z):

    M0: x* + x1,       z* + Dpsi zeta1
    M1: x* + x1 + x2,  z* + Dpsi (zeta1 + (x2, 0)) + 1/2 D^2 psi[zeta1, zeta1]
    M2: x* + x1 + x2', z* + Dpsi (zeta1 + (x2', 0))
        (x2' driven without KCL curvature)

The reference is the full nonlinear DAE (BDF, rtol 1e-10, atol 1e-12). The
observables are the COI frequency deviation, abs(V) at bus 16 and the current
magnitude of the converter at bus 30 (or of the machine at 30 when bus 30 is
not replaced). The error max_t abs(y_model − y_ref) is fitted as C eps^p over
eps in a geometric sequence. The discretization is exact for x1 (zero-order hold
on a piecewise-constant input) and first-order hold for the x2 forcing
(dt = 1e-3 s).

Expected from Taylor's theorem, and not claimed as new: p(M0) about 2 and
p(M1) about 3. The question is p(M2), and whether the KCL curvature is material.

Correction (v2, 2026-09-11). The first run (kept in v1_foh_bug/) had two defects,
found before any interpretation:

1. The first-order-hold coefficient was Gamma1 / dt^2 instead of Gamma1 / dt: the
   Van Loan block exponential already returns Gamma1 / dt. On a test ODE the error
   was 0.40 against 4.5e-7 after the fix. It corrupts x2 (M1, M2) only; x1 (M0)
   uses the zero-order hold and was correct.
2. At t = PULSE the reference evaluated the algebraic outputs with the pulse on
   and the models with it off. The algebraic outputs jump there, so that one
   sample is excluded from the comparison. The states are continuous.

v2 also parallelizes over (case, amplitude) and solves the reference's algebraic
equations by a chord Newton iteration (LU of G_z at the equilibrium, same stop
criterion max abs(g) < 1e-12, full Newton as fallback). Only speed changes.

v2 adds one control, M2t: M2 driven by the transverse part of x1 (component in
C = span{R_x, w} removed). With the full curvature this removal is exact
(D^2 r[zeta, c] = 0 for c in C, from the exact symmetries); the diagnostic
`symmetry_check_x2_full_vs_transverse` measures it. Without the KCL curvature
the identity fails, so M2 - M2t is the part of M2's error caused by the broken
rotation / drift symmetry.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _fc import WORKERS, FCExperiment, out_dir, write_json
from scipy.integrate import solve_ivp
from scipy.linalg import expm, lu_factor, lu_solve

from _f7_common import Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402

OUT = out_dir("FC07_second_order_curvature")
BUS = 16
PULSE = 0.5
T_END = 4.0
DT = 1e-3
EPS_MW = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]
CASES = [
    ("BASE P4", (), (0.03625, 1.425)),
    ("single 30 P4", (30,), (0.03625, 1.425)),
    ("critical triple 30+33+35 P4", (30, 33, 35), (0.03625, 1.425)),
    ("flagship P4 (transversely unstable)", (30, 33, 35, 37), (0.03625, 1.425)),
    ("flagship P_inf (safe policy)", (30, 33, 35, 37), (1.0, 0.5)),
]
_C: dict = {}


def _init():
    pin_blas_threads()
    _C["cfg"] = {c["name"]: c for c in r_configs()}["R_none"]


class Model:
    """F(x, z, u), G(x, z, u) with u = dP_load(bus 16) [pu], on the frozen DAE."""

    def __init__(self, case):
        self.case = case
        self.dae = case.dae
        self.net = self.dae.network
        self.p0 = self.net.loads[BUS]
        self.x0, self.z0 = case.equilibrium.x, case.equilibrium.z

    def _dae(self, u):
        from dataclasses import replace

        loads = dict(self.net.loads)
        loads[BUS] = self.p0 + complex(u, 0.0)
        return replace(self.dae, network=replace(self.net, loads=loads))

    def F(self, x, z, u):
        return self._dae(u).f(x, z, {}) if u else self.dae.f(x, z, {})

    def G(self, x, z, u):
        return self._dae(u).g(x, z, {}) if u else self.dae.g(x, z, {})


def outputs(model, x, z):
    dae = model.dae
    v = dae.voltages(z)
    om, m = [], []
    cur = None
    for slot in dae.slots:
        xs = x[slot.start : slot.stop]
        if slot.kind == "sg":
            om.append(xs[1])
            m.append(slot.device.parameters.m * slot.device.weight)
        if slot.bus == 30:
            vb = complex(v[dae.network.position(30)])
            cur = abs(slot.device.injection(xs, vb))
    coi = float(np.dot(m, om) / np.sum(m))
    return np.array([(coi - 1.0) * 60.0, abs(v[dae.network.position(BUS)]), cur])


def build_operators(model):
    x0, z0 = model.x0, model.z0
    jac = central_difference_jacobians(model.dae, x0, z0, {})
    lu = lu_factor(jac.gz)
    h = 1e-6
    fu = (model.F(x0, z0, h) - model.F(x0, z0, -h)) / (2 * h)
    gu = (model.G(x0, z0, h) - model.G(x0, z0, -h)) / (2 * h)
    jpsi = -lu_solve(lu, jac.gx)
    psiu = -lu_solve(lu, gu)
    a = jac.fx + jac.fz @ jpsi
    b = fu + jac.fz @ psiu
    # orthonormal basis of the exact center subspace C = span{R_x, w}
    r_x, _ = rotation_generator(model.dae, z0)
    w = frequency_partner(model.dae).w
    basis = np.column_stack([r_x] + ([] if w is None else [w]))
    uc = np.linalg.qr(basis)[0]
    return {
        "fz": jac.fz,
        "lu": lu,
        "jpsi": jpsi,
        "psiu": psiu,
        "a": a,
        "b": b,
        "uc": uc,
    }


def curvature(model, ops, xa, ua, hstep):
    """(D^2F[Xi zeta, Xi zeta], D^2 psi[zeta, zeta]) for zeta = (xa, ua)."""

    x0, z0 = model.x0, model.z0
    za = ops["jpsi"] @ xa + ops["psiu"] * ua
    scale = max(np.linalg.norm(xa), abs(ua), 1e-300)
    h = hstep / scale
    fp = model.F(x0 + h * xa, z0 + h * za, h * ua)
    fm = model.F(x0 - h * xa, z0 - h * za, -h * ua)
    gp = model.G(x0 + h * xa, z0 + h * za, h * ua)
    gm = model.G(x0 - h * xa, z0 - h * za, -h * ua)
    f0 = model.F(x0, z0, 0.0)
    g0 = model.G(x0, z0, 0.0)
    d2f = (fp - 2 * f0 + fm) / (h * h)
    d2g = (gp - 2 * g0 + gm) / (h * h)
    d2psi = -lu_solve(ops["lu"], d2g)
    return d2f, d2psi


def propagators(a, dt):
    n = a.shape[0]
    big = np.zeros((3 * n, 3 * n))
    big[:n, :n] = a * dt
    big[:n, n : 2 * n] = np.eye(n) * dt
    big[n : 2 * n, 2 * n :] = np.eye(n)
    e = expm(big)
    phi = e[:n, :n]
    g0 = e[:n, n : 2 * n]  # Gamma0 = int_0^dt e^{A s} ds
    g1 = e[:n, 2 * n :]  # Gamma1 / dt, Gamma1 = int_0^dt e^{A s} (dt - s) ds
    return phi, g0, g1


def solve_z(dae, x, z, lu):
    """Algebraic equations g(x, z) = 0 to max abs(g) < 1e-12: chord Newton with a
    fixed LU of G_z, then full Newton from the best iterate if the chord stalls."""

    best, best_r = z, np.inf
    for _ in range(60):
        r = dae.g(x, z, {})
        nr = np.abs(r).max()
        if nr < 1e-12:
            return z
        if nr < best_r:
            best, best_r = z, nr
        elif nr > 10 * best_r:
            break
        z = z - lu_solve(lu, r)
    z = best
    for _ in range(30):
        r = dae.g(x, z, {})
        if np.abs(r).max() < 1e-12:
            break
        j = central_difference_jacobians(dae, x, z, {})
        z = z - np.linalg.solve(j.gz, r)
    return z


def reference(model, eps_pu, times):
    dae_on = model._dae(eps_pu)
    dae_off = model.dae
    x0, z0 = model.x0.copy(), model.z0.copy()
    zstate = {"z": z0.copy()}
    lus = {
        id(d): lu_factor(central_difference_jacobians(d, x0, z0, {}).gz)
        for d in (dae_on, dae_off)
    }

    def rhs_for(dae):
        def rhs(t, x):
            z = solve_z(dae, x, zstate["z"], lus[id(dae)])
            zstate["z"] = z
            return dae.f(x, z, {})

        return rhs

    jac0 = central_difference_jacobians(model.dae, x0, z0, {})
    a0 = jac0.fx - jac0.fz @ np.linalg.solve(jac0.gz, jac0.gx)
    s1 = solve_ivp(
        rhs_for(dae_on),
        (0, PULSE),
        x0,
        method="BDF",
        rtol=1e-10,
        atol=1e-12,
        jac=lambda t, x: a0,
        dense_output=True,
        max_step=0.01,
    )
    s2 = solve_ivp(
        rhs_for(dae_off),
        (PULSE, T_END),
        s1.y[:, -1],
        method="BDF",
        rtol=1e-10,
        atol=1e-12,
        jac=lambda t, x: a0,
        dense_output=True,
        max_step=0.01,
    )
    ys = []
    for t in times:
        dae, sol = (dae_on, s1) if t <= PULSE else (dae_off, s2)
        x = sol.sol(t)
        z = solve_z(dae, x, zstate["z"], lus[id(dae)])
        zstate["z"] = z
        ys.append(outputs(model, x, z))
    return np.array(ys)


def case_task(task):
    (name, members, (g, k)), eps_mw = task
    case = solve_config(members, Theta(g, k, 1.5, 1.0), _C["cfg"])
    model = Model(case)
    ops = build_operators(model)
    a, b = ops["a"], ops["b"]
    n = a.shape[0]
    times = np.arange(0.0, T_END + DT / 2, DT)
    phi, g0, g1 = propagators(a, DT)
    obs_idx = np.arange(0, len(times), 20)
    # the algebraic outputs jump at the pulse end; compare the states elsewhere only
    obs_idx = obs_idx[np.abs(times[obs_idx] - PULSE) > DT / 2]
    eps = eps_mw / 100.0
    u = np.where(times < PULSE - 1e-12, eps, 0.0)
    x1 = np.zeros((len(times), n))
    for i in range(len(times) - 1):
        x1[i + 1] = phi @ x1[i] + g0 @ (b * u[i])
    # Second-order forcing along the trajectory. The exact symmetries give
    # D^2 r[zeta, c] = 0 for c in C, so the full forcing may equally be evaluated
    # on the transverse part of x1 (a consistency check). Without the KCL
    # curvature that identity fails; M2t shows how much of M2's error is the
    # broken symmetry acting on the common angle / frequency drift.
    nz = ops["fz"].shape[1]
    f = {k: np.zeros((len(times), n)) for k in ("full", "nokcl", "full_t", "nokcl_t")}
    d2psi_t = np.zeros((len(times), nz))
    kcl_ratio = []
    for i in range(len(times)):
        if not np.any(x1[i]) and u[i] == 0:
            continue
        d2f, d2psi = curvature(model, ops, x1[i], u[i], 1e-4)
        xt = x1[i] - ops["uc"] @ (ops["uc"].T @ x1[i])
        d2f_t, d2psi_tr = curvature(model, ops, xt, u[i], 1e-4)
        f["full"][i] = 0.5 * (d2f + ops["fz"] @ d2psi)
        f["nokcl"][i] = 0.5 * d2f
        f["full_t"][i] = 0.5 * (d2f_t + ops["fz"] @ d2psi_tr)
        f["nokcl_t"][i] = 0.5 * d2f_t
        d2psi_t[i] = d2psi
        kcl_ratio.append(
            np.linalg.norm(ops["fz"] @ d2psi) / max(np.linalg.norm(d2f), 1e-300)
        )
    x2 = {k: np.zeros((len(times), n)) for k in f}
    for k, fk in f.items():
        for i in range(len(times) - 1):
            x2[k][i + 1] = phi @ x2[k][i] + g0 @ fk[i] + g1 @ (fk[i + 1] - fk[i])
    ref = reference(model, eps, times[obs_idx])
    ys = {m: [] for m in ("M0", "M1", "M2", "M2t")}
    for i in obs_idx:
        z1 = ops["jpsi"] @ x1[i] + ops["psiu"] * u[i]
        xa, za = model.x0 + x1[i], model.z0 + z1
        ys["M0"].append(outputs(model, xa, za))
        ys["M1"].append(
            outputs(
                model,
                xa + x2["full"][i],
                za + ops["jpsi"] @ x2["full"][i] + 0.5 * d2psi_t[i],
            )
        )
        for m, k in (("M2", "nokcl"), ("M2t", "nokcl_t")):
            ys[m].append(outputs(model, xa + x2[k][i], za + ops["jpsi"] @ x2[k][i]))
    labels = {
        "M0": "M0 linear",
        "M1": "M1 full second order",
        "M2": "M2 no KCL curvature",
        "M2t": "M2t no KCL curvature, transverse forcing",
    }
    rows = []
    for m, y in ys.items():
        err = np.abs(np.array(y) - ref).max(axis=0)
        rows.append(
            {
                "case": name,
                "model": labels[m],
                "eps_mw": eps_mw,
                "err_freq_hz": err[0],
                "err_v16_pu": err[1],
                "err_i30_pu": err[2],
                "ref_peak_freq_hz": float(np.abs(ref[:, 0]).max()),
                "ref_peak_dv16_pu": float(np.abs(ref[:, 1] - ref[0, 1]).max()),
            }
        )
    scale = max(np.abs(x2["full"]).max(), 1e-300)
    return rows, {
        "case": name,
        "eps_mw": eps_mw,
        "kcl_to_f_ratio_median": float(np.median(kcl_ratio)),
        "kcl_to_f_ratio_max": float(np.max(kcl_ratio)),
        "symmetry_check_x2_full_vs_transverse": float(
            np.abs(x2["full"] - x2["full_t"]).max() / scale
        ),
        "x2_nokcl_over_full": float(np.abs(x2["nokcl"]).max() / scale),
        "x2_nokcl_t_over_full": float(np.abs(x2["nokcl_t"]).max() / scale),
    }


def slopes(frame):
    out = []
    for (case, model), blk in frame.groupby(["case", "model"]):
        rec = {"case": case, "model": model}
        for col in ("err_freq_hz", "err_v16_pu", "err_i30_pu"):
            e = blk[col].to_numpy()
            x = blk.eps_mw.to_numpy()
            ok = e > 0
            rec[f"p_{col}"] = (
                float(np.polyfit(np.log(x[ok]), np.log(e[ok]), 1)[0])
                if ok.sum() > 2
                else np.nan
            )
            # slope on the three smallest amplitudes (asymptotic regime)
            rec[f"p_small_{col}"] = float(
                np.polyfit(np.log(x[:4]), np.log(e[:4]), 1)[0]
            )
        out.append(rec)
    return pd.DataFrame(out)


def main(argv) -> int:
    exp = FCExperiment(
        name="FC07_second_order_curvature",
        question=(
            "Is the AC algebraic-manifold curvature needed for a consistent "
            "second-order response?"
        ),
        config={
            "bus": BUS,
            "pulse_s": PULSE,
            "t_end": T_END,
            "dt": DT,
            "eps_mw": EPS_MW,
            "cases": [c[0] for c in CASES],
        },
        workers=WORKERS,
    )
    started = time.time()
    tasks = [(c, e) for c in CASES for e in EPS_MW]
    tasks.sort(key=lambda t: -t[1])  # the large amplitudes (slowest references) first
    with Pool(min(WORKERS, len(tasks)), initializer=_init) as pool:
        res = pool.map(case_task, tasks, chunksize=1)
    frame = pd.DataFrame([r for rows, _ in res for r in rows])
    frame = frame.sort_values(["case", "model", "eps_mw"]).reset_index(drop=True)
    frame.to_csv(OUT / "FC07_errors.csv", index=False)
    sl = slopes(frame)
    sl.to_csv(OUT / "FC07_slopes.csv", index=False)
    ratios = [r for _, r in res]
    pd.DataFrame(ratios).to_csv(OUT / "FC07_diagnostics.csv", index=False)
    summary = {
        "version": "v2 (FOH coefficient and pulse-end sampling corrected)",
        "slopes": sl.to_dict("records"),
        "kcl_ratio": ratios,
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC07_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    pd.set_option("display.width", 250)
    print(sl.to_string(index=False))
    print(json.dumps(ratios, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
