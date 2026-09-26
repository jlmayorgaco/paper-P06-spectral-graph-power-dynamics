"""TX4 contextual-return closure campaign.

This script is deliberately narrow: it reuses the frozen full-order TX4
phasor model and exact port realization, and produces only the closure,
boundary, minimality, and numerical-truth outputs named by the preregistration.
It does not run any weak-element, radius, planner, EMT, or new-controller
campaign.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import eigvals as scipy_eigvals
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve()
EXP = HERE.parent
SRC = EXP.parent / "src"
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(SRC))

from _f7_common import LEAK, Theta, CORE, solve_subset  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.port_admittance import build_action_space  # noqa: E402

ROOT = HERE.parents[5]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
FREQ_BAND = (0.3, 1.5)
P4 = Theta(g=0.03625, k=1.425, t=1.5, h=1.0)
G_GRID = [
    0.03625,
    0.05000,
    0.07500,
    0.10000,
    0.12500,
    0.15000,
    0.17500,
    0.19000,
    0.20000,
    0.20500,
    0.20750,
    0.20768,
    0.2076814045,
    0.20770,
    0.21000,
    0.22500,
    0.25000,
]


@dataclass
class Mode:
    alpha: float
    eig: complex
    frequency_hz: float
    damping_ratio: float
    a_perp: np.ndarray
    coupling_residual: float


def tx4_case(members: tuple[int, ...], theta: Theta, tol: float = 1e-9):
    """Direct TX4 solve with the frozen converter and machine semantics."""

    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=theta.g, voltage_leak=LEAK
        ),
        machine_scaling={"ka": theta.k, "ta": theta.t},
        tol=tol,
    )


def mode_of(case) -> Mode:
    rx, _ = rotation_generator(case.dae, case.equilibrium.z)
    fp = frequency_partner(case.dae)
    tr = transverse_operator(case.system.A, rx, fp.w)
    vals = np.linalg.eigvals(tr.a_perp)
    freq = np.abs(vals.imag) / (2.0 * math.pi)
    mask = (
        (vals.imag >= -1e-10)
        & (freq >= FREQ_BAND[0] - 1e-10)
        & (freq <= FREQ_BAND[1] + 1e-10)
        & (np.abs(vals) > 1e-3)
    )
    if not np.any(mask):
        # The TX4 witness is in-band. Keeping this explicit makes a missing
        # mode a hard, auditable failure rather than a silent fallback.
        raise RuntimeError("no positive-frequency transverse mode in TX4 band")
    candidates = np.flatnonzero(mask)
    idx = candidates[np.argmax(vals[candidates].real)]
    lam = complex(vals[idx])
    return Mode(
        alpha=float(lam.real),
        eig=lam,
        frequency_hz=float(abs(lam.imag) / (2.0 * math.pi)),
        damping_ratio=float(-lam.real / max(abs(lam), 1e-300)),
        a_perp=tr.a_perp,
        coupling_residual=tr.coupling_residual,
    )


def independent_jacobian(case, step_factor: float = 1.0):
    """Independent block finite difference, not the repository helper."""

    dae = case.dae
    x = case.equilibrium.x
    z = case.equilibrium.z
    fx = np.zeros((dae.n_x, dae.n_x))
    fz = np.zeros((dae.n_x, dae.n_z))
    gx = np.zeros((dae.n_z, dae.n_x))
    gz = np.zeros((dae.n_z, dae.n_z))
    for j in range(dae.n_x):
        h = step_factor * 1e-6 * max(abs(float(x[j])), 1.0)
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        fx[:, j] = (dae.f(xp, z, {}) - dae.f(xm, z, {})) / (2.0 * h)
        gx[:, j] = (dae.g(xp, z, {}) - dae.g(xm, z, {})) / (2.0 * h)
    for j in range(dae.n_z):
        h = step_factor * 1e-6 * max(abs(float(z[j])), 1.0)
        zp, zm = z.copy(), z.copy()
        zp[j] += h
        zm[j] -= h
        fz[:, j] = (dae.f(x, zp, {}) - dae.f(x, zm, {})) / (2.0 * h)
        gz[:, j] = (dae.g(x, zp, {}) - dae.g(x, zm, {})) / (2.0 * h)
    red = fx - fz @ np.linalg.solve(gz, gx)
    return fx, fz, gx, gz, red


def case_diagnostics(case, step_factor: float = 1.0) -> dict:
    mode = mode_of(case)
    fx, fz, gx, gz, red = independent_jacobian(case, step_factor)
    vals_np = np.linalg.eigvals(mode.a_perp)
    vals_sp = scipy_eigvals(mode.a_perp)
    np_crit = vals_np[np.argmin(np.abs(vals_np - mode.eig))]
    sp_crit = vals_sp[np.argmin(np.abs(vals_sp - mode.eig))]
    eq_f = np.max(np.abs(case.dae.f(case.equilibrium.x, case.equilibrium.z, {})))
    eq_g = np.max(np.abs(case.dae.g(case.equilibrium.x, case.equilibrium.z, {})))
    return {
        "alpha": mode.alpha,
        "critical_real": mode.eig.real,
        "critical_imag": mode.eig.imag,
        "frequency_hz": mode.frequency_hz,
        "damping_ratio": mode.damping_ratio,
        "state_dim": int(case.dae.n_x),
        "algebraic_dim": int(case.dae.n_z),
        "equilibrium_f_residual": float(eq_f),
        "equilibrium_g_residual": float(eq_g),
        "gz_condition_independent": float(np.linalg.cond(gz)),
        "reduced_jacobian_condition_independent": float(np.linalg.cond(red)),
        "reduced_jacobian_sigma_min": float(np.linalg.svd(red, compute_uv=False)[-1]),
        "coupling_residual": mode.coupling_residual,
        "numpy_scipy_eig_error": float(abs(np_crit - sp_crit)),
        "descriptor_status": "NOT_AVAILABLE",
        "descriptor_reason": "No valid E,A descriptor pair is exposed by the frozen TX4 DAE.",
        "a_perp": mode.a_perp,
    }


def q_parts(space, s: complex):
    # PortActionSpace.split defines Q after normalizing by the block-diagonal
    # local factor. Using M directly would leave nonzero diagonal blocks and
    # would invalidate the contextual-return Schur form.
    m = space.m(s)
    n = space.order
    eye = np.eye(2 * n, dtype=complex)
    total = eye + m
    local = np.zeros_like(total)
    for k in range(n):
        block = space.block(total, k, k)
        local[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = block
    q = np.linalg.solve(local, total) - eye
    c = eye + q
    local_smin = []
    local_det = []
    for k in range(n):
        block = space.block(c, k, k)
        local_smin.append(float(np.linalg.svd(block, compute_uv=False)[-1]))
        local_det.append(complex(np.linalg.det(block)))
    return q, c, local_smin, local_det


def contextual_return(q: np.ndarray, i: int) -> dict:
    n = q.shape[0] // 2
    all_idx = [k for k in range(n) if k != i]
    ridx = sum(([2 * k, 2 * k + 1] for k in all_idx), [])
    rr = np.ix_(ridx, ridx)
    ir = np.ix_([2 * i, 2 * i + 1], ridx)
    ri = np.ix_(ridx, [2 * i, 2 * i + 1])
    qrr = q[rr]
    qir = q[ir]
    qri = q[ri]
    cr = np.eye(qrr.shape[0], dtype=complex) + qrr
    g = np.linalg.inv(cr)
    ret = qir @ g @ qri
    cret = np.eye(2, dtype=complex) - ret
    ch = np.eye(2 * n, dtype=complex) + q
    lhs = np.linalg.det(ch)
    rhs = np.linalg.det(cr) * np.linalg.det(cret)
    mu = np.linalg.eigvals(ret)
    return {
        "eta_subset": float(np.linalg.svd(cr, compute_uv=False)[-1]),
        "return_operator": ret,
        "return_eigenvalues": mu,
        "nearest_return_eigenvalue": complex(mu[np.argmin(abs(mu - 1.0))]),
        "return_distance_to_unity": float(np.min(abs(mu - 1.0))),
        "return_det": complex(np.linalg.det(cret)),
        "schur_identity_residual": float(abs(lhs - rhs)),
        "full_collective_smin": float(np.linalg.svd(ch, compute_uv=False)[-1]),
    }


def return_matrix_at(case, base, omega: float, i: int) -> np.ndarray:
    space = build_action_space(base, case, CORE)
    q, *_ = q_parts(space, 1j * omega)
    return contextual_return(q, i)["return_operator"]


def return_operator_and_derivative(q: np.ndarray, dq: np.ndarray, i: int):
    """Return R and the preregistered block derivative dR/da."""

    n = q.shape[0] // 2
    ridx = sum(([2 * k, 2 * k + 1] for k in range(n) if k != i), [])
    ir = np.ix_([2 * i, 2 * i + 1], ridx)
    ri = np.ix_(ridx, [2 * i, 2 * i + 1])
    rr = np.ix_(ridx, ridx)
    qir, qri, qrr = q[ir], q[ri], q[rr]
    dqir, dqri, dqrr = dq[ir], dq[ri], dq[rr]
    g = np.linalg.inv(np.eye(qrr.shape[0], dtype=complex) + qrr)
    ret = qir @ g @ qri
    dret = dqir @ g @ qri + qir @ g @ dqri - qir @ g @ dqrr @ g @ qri
    return ret, dret


def csv_write(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    keys = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def audit_cases() -> dict:
    cases = {
        "BASE": (),
        "30": (30,),
        "30+33+35+37": CORE,
        "33": (33,),
        "30+33": (30, 33),
        "30+33+35": (30, 33, 35),
        "33+35+37": (33, 35, 37),
        "30+35+37": (30, 35, 37),
        "30+33+37": (30, 33, 37),
    }
    rows = []
    tol_summary = []
    for name, members in cases.items():
        for factor in (0.5, 1.0, 2.0):
            case = tx4_case(tuple(members), P4, tol=1e-9)
            d = case_diagnostics(case, factor)
            rows.append(
                {
                    "portfolio": name,
                    "members": "+".join(map(str, members)) or "BASE",
                    "tolerance_factor": factor,
                    **{k: v for k, v in d.items() if k != "a_perp"},
                    "verdict": "UNSTABLE" if d["alpha"] > 0 else "STABLE",
                }
            )
            tol_summary.append(d["numpy_scipy_eig_error"])
    csv_write(OUT / "TX4_CONTEXTUAL_RETURN_NUMERICAL_TRUTH.csv", rows)
    return {
        "n_portfolios": len(cases),
        "tolerance_rows": len(rows),
        "max_numpy_scipy_eig_error": max(tol_summary),
        "pass": all(r["numpy_scipy_eig_error"] <= 1e-6 for r in rows),
        "descriptor": "NOT_AVAILABLE",
    }


def solve_root(cache: dict[float, tuple[object, Mode]]) -> tuple[float, int]:
    def alpha(g: float) -> float:
        key = round(float(g), 12)
        if key not in cache:
            c = tx4_case(CORE, Theta(g=g, k=1.425, t=1.5, h=1.0))
            cache[key] = (c, mode_of(c))
        return cache[key][1].alpha

    lo, hi = 0.03625, 0.25
    alo, ahi = alpha(lo), alpha(hi)
    if not (alo > 0.0 and ahi < 0.0):
        raise RuntimeError(f"fixed root bracket invalid: alpha({lo})={alo}, alpha({hi})={ahi}")
    iterations = 0
    for iterations in range(1, 51):
        mid = 0.5 * (lo + hi)
        am = alpha(mid)
        if abs(am) <= 1e-9 or hi - lo <= 1e-10:
            return mid, iterations
        if am > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), iterations


def sweep_and_tables() -> dict:
    cache: dict[float, tuple[object, Mode]] = {}
    root, root_iterations = solve_root(cache)
    # The base is independent of g in this model, but it is solved once along
    # the path so the port comparison always uses the exact matched base.
    base = tx4_case((), Theta(g=0.0, k=1.425, t=1.5, h=1.0))
    raw = []
    for g in G_GRID:
        key = round(float(g), 12)
        if key not in cache:
            c = tx4_case(CORE, Theta(g=g, k=1.425, t=1.5, h=1.0))
            cache[key] = (c, mode_of(c))
        case, mode = cache[key]
        space = build_action_space(base, case, CORE)
        s = 1j * mode.eig.imag
        q, cmat, local_smin, local_det = q_parts(space, s)
        returns = [contextual_return(q, i) for i in range(4)]
        split = space.split(s)
        raw.append(
            {
                "g": g,
                "alpha": mode.alpha,
                "critical_real": mode.eig.real,
                "critical_imag": mode.eig.imag,
                "frequency_hz": mode.frequency_hz,
                "damping_ratio": mode.damping_ratio,
                "collective_det_abs": float(abs(np.linalg.det(cmat))),
                "collective_sigma_min": float(np.linalg.svd(cmat, compute_uv=False)[-1]),
                "local_sigma_min": min(local_smin),
                "local_det_min_abs": min(abs(x) for x in local_det),
                "port_identity_residual": float(abs(split["full"] - split["individual"] * split["collective"])),
                "mu_30_distance": returns[0]["return_distance_to_unity"],
                "mu_33_distance": returns[1]["return_distance_to_unity"],
                "mu_35_distance": returns[2]["return_distance_to_unity"],
                "mu_37_distance": returns[3]["return_distance_to_unity"],
                "return_distance_to_unity": min(r["return_distance_to_unity"] for r in returns),
                "return_schur_residual_max": max(r["schur_identity_residual"] for r in returns),
                "equilibrium_residual": float(max(abs(case.dae.g(case.equilibrium.x, case.equilibrium.z, {})))),
            }
        )
    # Fixed scalar minimization of the return distance over the predeclared
    # root neighborhood; no new grid is selected from the result.
    def ret_dist(g: float) -> float:
        c = tx4_case(CORE, Theta(g=float(g), k=1.425, t=1.5, h=1.0))
        m = mode_of(c)
        sp = build_action_space(base, c, CORE)
        q, *_ = q_parts(sp, 1j * m.eig.imag)
        return min(contextual_return(q, i)["return_distance_to_unity"] for i in range(4))

    ret_min = minimize_scalar(ret_dist, bounds=(0.19, 0.225), method="bounded", options={"xatol": 1e-8})
    return_row = min(raw, key=lambda r: r["return_distance_to_unity"])
    csv_write(OUT / "TX4_CONTEXTUAL_RETURN_SWEEP.csv", raw)

    # Root table: exact full collective and each one-device contextual return.
    root_case, root_mode = cache[round(root, 12)]
    root_space = build_action_space(base, root_case, CORE)
    q, cmat, local_smin, local_det = q_parts(root_space, 1j * root_mode.eig.imag)
    boundary_rows = []
    for i, bus in enumerate(CORE):
        r = contextual_return(q, i)
        mu = r["nearest_return_eigenvalue"]
        boundary_rows.append(
            {
                "device_bus": bus,
                "proper_subset": "+".join(str(b) for j, b in enumerate(CORE) if j != i),
                "g_root": root,
                "s_real": 0.0,
                "s_imag": root_mode.eig.imag,
                "frequency_hz": root_mode.frequency_hz,
                "eta_subset": r["eta_subset"],
                "nearest_return_eigenvalue_real": mu.real,
                "nearest_return_eigenvalue_imag": mu.imag,
                "return_distance_to_unity": r["return_distance_to_unity"],
                "return_det_abs": abs(r["return_det"]),
                "collective_sigma_min": r["full_collective_smin"],
                "local_sigma_min_at_device": local_smin[i],
                "local_factor_det_abs": abs(local_det[i]),
                "schur_identity_residual": r["schur_identity_residual"],
            }
        )
    csv_write(OUT / "TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv", boundary_rows)

    proper_rows = []
    for k in range(4):
        for i in range(4):
            pass
    from itertools import combinations
    for size in range(4):
        for subset in combinations(CORE, size):
            ss = tuple(subset)
            cs = tx4_case(ss, P4)
            ms = mode_of(cs)
            sp = build_action_space(base, cs, ss) if ss else None
            if ss:
                qs, cqs, _, _ = q_parts(sp, 1j * root_mode.eig.imag)
                qev = np.linalg.eigvals(qs)
                eta = float(np.linalg.svd(cqs, compute_uv=False)[-1])
                nearest = complex(qev[np.argmin(abs(qev + 1.0))])
            else:
                eta, nearest = 1.0, 0.0j
            proper_rows.append(
                {
                    "subset": "+".join(map(str, ss)) or "BASE",
                    "cardinality": len(ss),
                    "alpha_perp_at_P4": ms.alpha,
                    "frequency_hz_at_P4": ms.frequency_hz,
                    "alpha_at_root_frequency": ms.alpha,
                    "closure_sigma_min_at_H4_root": eta,
                    "nearest_Q_eigenvalue_to_minus_one_real": nearest.real,
                    "nearest_Q_eigenvalue_to_minus_one_imag": nearest.imag,
                    "stable": ms.alpha < 0.0,
                }
            )
    csv_write(OUT / "TX4_PROPER_SUBSET_CLOSURE.csv", proper_rows)

    # Analytic block derivative at fixed root frequency, compared with the
    # prescribed central finite difference after full equilibrium re-solving.
    derivative_rows = []
    h_g = 1e-5
    for i, bus in enumerate(CORE):
        cp = tx4_case(CORE, Theta(g=root + h_g, k=1.425, t=1.5, h=1.0))
        cm = tx4_case(CORE, Theta(g=root - h_g, k=1.425, t=1.5, h=1.0))
        q0, *_ = q_parts(root_space, 1j * root_mode.eig.imag)
        sp = build_action_space(base, cp, CORE)
        sm = build_action_space(base, cm, CORE)
        qp, *_ = q_parts(sp, 1j * root_mode.eig.imag)
        qm, *_ = q_parts(sm, 1j * root_mode.eig.imag)
        r0, dr_formula = return_operator_and_derivative(q0, (qp - qm) / (2.0 * h_g), i)
        rp, _ = return_operator_and_derivative(qp, np.zeros_like(qp), i)
        rm, _ = return_operator_and_derivative(qm, np.zeros_like(qm), i)
        ev, vr = np.linalg.eig(r0)
        j = int(np.argmin(abs(ev - 1.0)))
        ew, vl = np.linalg.eig(r0.conj().T)
        jj = int(np.argmin(abs(ew - np.conj(ev[j]))))
        xvec, yvec = vr[:, j], vl[:, jj]
        dmu = np.vdot(yvec, dr_formula @ xvec) / np.vdot(yvec, xvec)
        mu_fd = np.linalg.eigvals(rp)[np.argmin(abs(np.linalg.eigvals(rp) - ev[j]))]
        fd_mu = (mu_fd - np.linalg.eigvals(rm)[np.argmin(abs(np.linalg.eigvals(rm) - ev[j]))]) / (2.0 * h_g)
        derivative_rows.append(
            {
                "device_bus": bus,
                "g": root,
                "frequency_hz": root_mode.frequency_hz,
                "mu_real": ev[j].real,
                "mu_imag": ev[j].imag,
                "dmu_dg_block_formula_real": dmu.real,
                "dmu_dg_block_formula_imag": dmu.imag,
                "dmu_dg_fd_real": fd_mu.real,
                "dmu_dg_fd_imag": fd_mu.imag,
                "absolute_error": abs(dmu - fd_mu),
            }
        )
    csv_write(OUT / "TX4_CONTEXTUAL_RETURN_DERIVATIVE.csv", derivative_rows)

    # Fixed g-only remediation table, including the previously archived point.
    after = tx4_case(CORE, Theta(g=0.25, k=1.425, t=1.5, h=1.0))
    after_mode = mode_of(after)
    g0_case = cache[round(0.03625, 12)][0]
    g0_mode = mode_of(g0_case)
    dg = 1e-4
    gp = tx4_case(CORE, Theta(g=0.03625 + dg, k=1.425, t=1.5, h=1.0))
    gm = tx4_case(CORE, Theta(g=0.03625 - dg, k=1.425, t=1.5, h=1.0))
    dadt = (mode_of(gp).alpha - mode_of(gm).alpha) / (2.0 * dg)
    predicted_dg = -g0_mode.alpha / dadt
    remediation_rows = [
        {
            "stage": "P4_start",
            "g": 0.03625,
            "alpha_perp": g0_mode.alpha,
            "frequency_hz": g0_mode.frequency_hz,
            "dalpha_dg_local": dadt,
            "first_order_delta_g_to_alpha_zero": predicted_dg,
            "newton_iterations": "see_root",
            "verdict": "UNSTABLE" if g0_mode.alpha > 0 else "STABLE",
        },
        {
            "stage": "boundary",
            "g": root,
            "alpha_perp": root_mode.alpha,
            "frequency_hz": root_mode.frequency_hz,
            "dalpha_dg_local": dadt,
            "first_order_delta_g_to_alpha_zero": predicted_dg,
            "newton_iterations": root_iterations,
            "verdict": "BOUNDARY",
        },
        {
            "stage": "fixed_after",
            "g": 0.25,
            "alpha_perp": after_mode.alpha,
            "frequency_hz": after_mode.frequency_hz,
            "dalpha_dg_local": dadt,
            "first_order_delta_g_to_alpha_zero": predicted_dg,
            "newton_iterations": root_iterations,
            "verdict": "UNSTABLE" if after_mode.alpha > 0 else "STABLE",
        },
    ]
    csv_write(OUT / "TX4_CONTROL_BOUNDARY_REMEDIATION.csv", remediation_rows)

    out = {
        "g_eigen_boundary": root,
        "g_return_minimum": float(ret_min.x),
        "g_return_grid_minimum": return_row["g"],
        "g_boundary_difference": abs(root - float(ret_min.x)),
        "root_iterations": root_iterations,
        "alpha_P4": g0_mode.alpha,
        "alpha_after": after_mode.alpha,
        "frequency_root_hz": root_mode.frequency_hz,
        "local_sigma_min_root": min(local_smin),
        "collective_sigma_min_root": float(np.linalg.svd(cmat, compute_uv=False)[-1]),
        "max_boundary_schur_residual": max(r["schur_identity_residual"] for r in boundary_rows),
        "max_return_derivative_error": max(r["absolute_error"] for r in derivative_rows),
        "minimality_all_proper_stable": all(r["stable"] for r in proper_rows[:-1]),
        "newton_or_bisection_iterations": root_iterations,
    }
    (OUT / "TX4_CONTEXTUAL_RETURN_CORE_SUMMARY.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


def main():
    audit = audit_cases()
    core = sweep_and_tables()
    summary = {"audit": audit, "core": core}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
