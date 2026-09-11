# ruff: noqa: E501  -- table labels kept on one line
"""PCV04 - Phase 6 remediation and Task E intervention ranking at P4 (complete H4).

Preregistration 5d0b1986, sections 5 (Task E) and 8. Documented intervention families
only: Q/V gain g, AVR gain scale k, damped condenser (frozen), documented governors
(frozen), the 12 frozen holdout lines. No common cost; |chi| is not an objective.

Port derivative: ds*/da = -(dD/da)/(dD/ds) with D(s, a) = det(I + M_SS(s; a)) on the
C2 realization, at the P4 critical eigenvalue of H4 (a simple zero of D).
Conventional sensitivity: central difference of the direct transverse eigenvalue.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _pcv import (
    CORE,
    FC03_POINTS,
    FC18_SUMMARY,
    H4,
    HOLDOUT,
    POINTS,
    ROOT,
    SCALE,
    Realization,
    kwargs_of,
    label,
    out_dir,
    transverse_of,
    write_json,
)
from scipy.stats import kendalltau, spearmanr

from _f7_common import Theta, solve_subset
from _overnight import pin_blas_threads
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator
from ibr_cycles.certification.transverse import transverse_operator
from ibr_cycles.diagnosis.baselines import generalized_scr
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

OUT = out_dir("PCV04")
PCV02 = OUT.parent / "PCV02"
P4 = POINTS["P4"]
NAMES = ("g", "k", "t", "h")
BOX = {"g": (0.0, 1.0), "k": (0.5, 2.3), "t": (0.5, 3.0), "h": (0.0, 2.0)}
STEP_FD = 1e-4
EPS_LINE = 0.002
GAMMAS = (1.25, 1.5, 2.0)
F2_STATIC = (
    ROOT
    / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff"
    / "dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation/data/F2_static_vs_dynamic_line_metrics.csv"
)
payload = json.loads((ROOT / "configs/ias2026/ieee39_network.json").read_text())
NET = load_network(ROOT / "configs/ias2026/ieee39_network.json")


# ---- verbatim from dynamic_forest_line_port_holdout.py (frozen F2c) ----
def ybus_scaled(line_idx, gamma):
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {b: i for i, b in enumerate(order)}
    n = len(order)
    y = np.zeros((n, n), complex)
    for ell, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        scale = gamma if ell == line_idx else 1.0
        series = scale / complex(line["r"], line["x"])
        charging = scale * complex(line["g"], line["b"]) / 2
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        m2 = abs(m) ** 2
        y[f, f] += (series + charging) / m2
        y[t, t] += series + charging
        y[f, t] += -series / np.conj(m)
        y[t, f] += -series / m
    for sh in payload["shunts"]:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return y


# ---- end verbatim ----


def scaled_net(li, gamma):
    return replace(NET, ybus=ybus_scaled(li, gamma))


def transverse_eigs(case):
    from ibr_cycles.certification.physical import physical_matrices

    a, _, _ = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(a, r_x, frequency_partner(case.dae).w)
    return np.linalg.eigvals(tr.a_perp)


def d_value(r: Realization, s):
    d, k = r.m_matrix(s)
    return complex(np.linalg.det(np.eye(d.shape[0]) + d @ k))


# ------------------------------------------------------------------ workers --
def w_alpha(task):
    """Exact alpha_perp(H4) at a policy point, optionally with a scaled line."""

    pin_blas_threads()
    tag, theta, line = task
    if line is None:
        case = solve_subset(H4, Theta(*theta))
    else:
        li, gamma = line
        case = solve_case(
            ReplacementPlan.of({b: 1.0 for b in H4}),
            network=scaled_net(li, gamma),
            **kwargs_of(theta),
        )
    return tag, transverse_of(case)


def w_eigs(task):
    pin_blas_threads()
    tag, theta, line = task
    if line is None:
        case = solve_subset(H4, Theta(*theta))
    else:
        li, gamma = line
        case = solve_case(
            ReplacementPlan.of({b: 1.0 for b in H4}),
            network=scaled_net(li, gamma),
            **kwargs_of(theta),
        )
    return tag, transverse_eigs(case)


def w_dvals(task):
    """D(s*) on a perturbed C2 realization (policy or line), for the port derivative."""

    pin_blas_threads()
    tag, theta, line, s_star = task
    extra = {} if line is None else {"network": scaled_net(*line)}
    r = Realization(CORE, theta, **extra)
    return tag, d_value(r, s_star)


def w_center(task):
    pin_blas_threads()
    theta, lam = task
    r = Realization(CORE, theta)
    ev = r.vertex_spectrum(H4)
    s = complex(ev[np.argmin(np.abs(ev - lam))])
    return s, d_value(r, s), d_value(r, s + 1e-5), d_value(r, s - 1e-5)


def port_derivative(theta, lam, coords):
    """Port ds*/da at the exact eigenvalue near lam for the named policy coordinates."""

    s, _, dp, dm = w_center((theta, lam))
    ds = (dp - dm) / 2e-5
    out = {}
    for name in coords:
        i = NAMES.index(name)
        tp, tm = list(theta), list(theta)
        tp[i] += STEP_FD
        tm[i] -= STEP_FD
        _, vp = w_dvals(("p", tuple(tp), None, s))
        _, vm = w_dvals(("m", tuple(tm), None, s))
        out[name] = complex(-((vp - vm) / (2 * STEP_FD)) / ds)
    return s, out


def critical_direct(theta):
    ev = transverse_eigs(solve_subset(H4, Theta(*theta)))
    f = ev.imag / (2 * np.pi)
    band = ev[(f >= 0.3) & (f <= 1.5)]
    return complex(band[np.argmax(band.real)])


def newton(coord, theta0, max_iter=8, tol=1e-7):
    pin_blas_threads()
    th = list(theta0)
    i = NAMES.index(coord)
    it = []
    for n in range(max_iter + 1):
        lam = critical_direct(tuple(th))
        s, grad = port_derivative(tuple(th), lam, (coord,))
        slope = grad[coord].real
        it.append(
            {
                "iter": n,
                coord: th[i],
                "re_s": s.real,
                "freq_hz": s.imag / (2 * np.pi),
                "dRe_s_d" + coord: slope,
            }
        )
        if abs(s.real) < tol:
            break
        th[i] += -s.real / slope
    return it


def bisect(coord, lo, hi, theta0, iters=40):
    pin_blas_threads()
    i = NAMES.index(coord)

    def a(v):
        th = list(theta0)
        th[i] = v
        return transverse_of(solve_subset(H4, Theta(*th)), classify=False)["alpha"]

    a_lo, a_hi = a(lo), a(hi)
    assert np.sign(a_lo) != np.sign(a_hi), (a_lo, a_hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if np.sign(a(mid)) == np.sign(a_lo):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def w_newton(task):
    return task[0], newton(*task)


def w_bisect(task):
    coord, lo, hi = task
    return coord, bisect(coord, lo, hi, P4)


# --------------------------------------------------------------------- main --
def main(argv) -> int:
    truth = pd.read_csv(PCV02 / "PCV02_truth_core.csv")
    t4 = truth[(truth.point == "P4") & (truth.subset == label(H4))].iloc[0]
    lam_p4 = complex(t4.crit_re, 2 * np.pi * t4.crit_hz)
    alpha_p4 = float(t4.alpha)
    holdout = json.loads(HOLDOUT.read_text())
    lines = [int(x) for x in holdout["lines"]]

    # ---------------- task lists ---------------------------------------------------
    alpha_tasks, eig_tasks, dval_tasks = [], [], []
    for name in NAMES:
        i = NAMES.index(name)
        for sgn in (+1, -1):
            v = P4[i] + sgn * 0.10 * SCALE[i]
            if BOX[name][0] <= v <= BOX[name][1]:
                th = list(P4)
                th[i] = v
                alpha_tasks.append((f"policy|{name}|{sgn}", tuple(th), None))
    for li in lines:
        for gm in GAMMAS:
            alpha_tasks.append((f"line|{li}|{gm}", P4, (li, gm)))
        for sgn, gm in ((+1, 1 + EPS_LINE), (-1, 1 - EPS_LINE)):
            eig_tasks.append((f"line|{li}|{sgn}", P4, (li, gm)))
    with Pool(16) as pool:
        center = pool.apply_async(w_center, ((P4, lam_p4),))
        alpha_res = dict(pool.map(w_alpha, alpha_tasks, chunksize=1))
        eig_res = dict(pool.map(w_eigs, eig_tasks, chunksize=1))
        s_star, _, dp, dm = center.get()
        ds = (dp - dm) / 2e-5
        for name in NAMES:
            i = NAMES.index(name)
            for sgn in (+1, -1):
                th = list(P4)
                th[i] += sgn * STEP_FD
                dval_tasks.append((f"policy|{name}|{sgn}", tuple(th), None, s_star))
        for li in lines:
            for sgn, gm in ((+1, 1 + EPS_LINE), (-1, 1 - EPS_LINE)):
                dval_tasks.append((f"line|{li}|{sgn}", P4, (li, gm), s_star))
        dval_res = dict(pool.map(w_dvals, dval_tasks, chunksize=1))
        newton_res = dict(pool.map(w_newton, [("g", P4), ("k", P4)], chunksize=1))
        bis_res = dict(
            pool.map(w_bisect, [("k", 1.30, 1.31), ("g", 0.2025, 0.2256)], chunksize=1)
        )

    # ---------------- policy coordinates (Task E) -----------------------------------
    prow = []
    for name in NAMES:
        i = NAMES.index(name)
        port = complex(
            -(
                (dval_res[f"policy|{name}|1"] - dval_res[f"policy|{name}|-1"])
                / (2 * STEP_FD)
            )
            / ds
        )
        best = None
        for sgn in (+1, -1):
            key = f"policy|{name}|{sgn}"
            if key in alpha_res:
                dA = alpha_res[key]["alpha"] - alpha_p4
                if best is None or dA < best[1]:
                    best = (sgn, dA, alpha_res[key]["alpha"], alpha_res[key]["status"])
        step = 0.10 * SCALE[i]
        prow.append(
            {
                "coordinate": name,
                "port_dRe_s": port.real,
                "port_dIm_s_hz": port.imag / (2 * np.pi),
                "step": step,
                "stabilizing_direction": best[0],
                "actual_delta_alpha": best[1],
                "alpha_after": best[2],
                "status_after": best[3],
                "predicted_delta_alpha_linear": -abs(port.real) * step,
                "predicted_direction": -int(np.sign(port.real)),
                "direction_correct": int(np.sign(port.real)) == -best[0],
            }
        )
    policy = pd.DataFrame(prow)
    policy.to_csv(OUT / "PCV04_policy_task_e.csv", index=False)

    # ---------------- lines (Task E + remediation) -----------------------------------
    f2 = pd.read_csv(F2_STATIC).set_index("line_index")
    flow0 = solve_power_flow(NET)
    gscr0 = generalized_scr(NET, H4, flow0)
    lrow = []
    for li in lines:
        port = complex(
            -((dval_res[f"line|{li}|1"] - dval_res[f"line|{li}|-1"]) / (2 * EPS_LINE))
            / ds
        )
        ep, em = eig_res[f"line|{li}|1"], eig_res[f"line|{li}|-1"]
        lp = ep[np.argmin(np.abs(ep - lam_p4))]
        lm = em[np.argmin(np.abs(em - lam_p4))]
        conv = complex((lp - lm) / (2 * EPS_LINE))
        net15 = scaled_net(li, 1.5)
        dgscr = generalized_scr(net15, H4, solve_power_flow(net15)) - gscr0
        rec = {
            "line_index": li,
            "line": f"L{li:02d}:{int(payload['lines'][li]['bus1'])}-{int(payload['lines'][li]['bus2'])}",
            "port_dRe_s_dgamma": port.real,
            "conventional_dRe_lambda_dgamma": conv.real,
            "pred_port_delta_alpha_1.5": 0.5 * port.real,
            "pred_conv_delta_alpha_1.5": 0.5 * conv.real,
            "delta_gscr_1.5": dgscr,
            "inverse_x": float(f2.loc[li, "inverse_x"]),
            "x_weighted_betweenness": float(f2.loc[li, "x_weighted_betweenness"]),
        }
        for gm in GAMMAS:
            a = alpha_res[f"line|{li}|{gm}"]
            rec[f"alpha_{gm}"] = a["alpha"]
            rec[f"status_{gm}"] = a["status"]
            rec[f"actual_delta_alpha_{gm}"] = a["alpha"] - alpha_p4
        lrow.append(rec)
    lines_df = pd.DataFrame(lrow)
    lines_df.to_csv(OUT / "PCV04_lines.csv", index=False)
    actual_stab = -lines_df["actual_delta_alpha_1.5"]
    te = []
    for name, pred in (
        ("port derivative (framework)", -lines_df["pred_port_delta_alpha_1.5"]),
        (
            "conventional full-DAE eigen-sensitivity (B3-type)",
            -lines_df["pred_conv_delta_alpha_1.5"],
        ),
        ("static: delta gSCR(H4) at 1.5x", lines_df["delta_gscr_1.5"]),
        ("static: 1/x", lines_df["inverse_x"]),
        (
            "static: x-weighted betweenness (frozen F2)",
            lines_df["x_weighted_betweenness"],
        ),
    ):
        te.append(
            {
                "intervention_set": "12 frozen lines at P4, gamma = 1.5",
                "method": name,
                "spearman": float(spearmanr(pred, actual_stab).statistic),
                "kendall": float(kendalltau(pred, actual_stab).statistic),
                "top1_agrees": int(np.argmax(pred.to_numpy()))
                == int(np.argmax(actual_stab.to_numpy())),
                "finite_step_sign_agreement": float(
                    (np.sign(pred) == np.sign(actual_stab)).mean()
                )
                if not name.startswith("static")
                else float("nan"),
            }
        )
    pa = -policy.actual_delta_alpha
    pp = -policy.predicted_delta_alpha_linear
    te.append(
        {
            "intervention_set": "policy coordinates g,k,t,h at P4 (n = 4)",
            "method": "port derivative (framework)",
            "spearman": float(spearmanr(pp, pa).statistic),
            "kendall": float(kendalltau(pp, pa).statistic),
            "top1_agrees": int(np.argmax(pp.to_numpy()))
            == int(np.argmax(pa.to_numpy())),
        }
    )
    task_e = pd.DataFrame(te)
    task_e.to_csv(OUT / "PCV04_taskE.csv", index=False)

    # ---------------- remediation table ------------------------------------------------
    fc18 = json.loads(FC18_SUMMARY.read_text())
    fc03 = pd.read_csv(FC03_POINTS)
    cond = fc03[fc03.tag == "CONDENSER P4"][
        ["condenser", "alpha_flag_frozen", "H_frozen", "alpha_flag_governed"]
    ]
    gov = fc03[fc03.tag == "F8 P4"].iloc[0]
    g_newton, k_newton = newton_res["g"], newton_res["k"]
    rem = [
        {
            "intervention": "Q/V gain g (converter policy)",
            "initial_alpha": alpha_p4,
            "final_state": f"g* = {bis_res['g']:.8f} (direct bisection); stable side G_S alpha from PCV02",
            "port_local_derivative": policy.set_index("coordinate").loc[
                "g", "port_dRe_s"
            ],
            "one_step_prediction": f"Delta g = {-alpha_p4 / policy.set_index('coordinate').loc['g', 'port_dRe_s']:.4f} vs actual {bis_res['g'] - P4[0]:.4f}",
            "newton": f"{len(g_newton) - 1} steps to g = {g_newton[-1]['g']:.10f} (|Re s| = {abs(g_newton[-1]['re_s']):.1e})",
            "frozen_reference": f"FC18 continuation g* = {fc18['retune_direct_boundary_g']}",
            "tds": "none at G_S (frozen G2 covers P4 unstable and P_inf stable)",
            "independent_tool": "no (custom GFL)",
        },
        {
            "intervention": "AVR gain scale k (SG excitation policy)",
            "initial_alpha": alpha_p4,
            "final_state": f"k* = {bis_res['k']:.8f} (direct bisection in [1.30, 1.31]); frozen F7A bracket [1.303125, 1.30625]",
            "port_local_derivative": policy.set_index("coordinate").loc[
                "k", "port_dRe_s"
            ],
            "one_step_prediction": f"Delta k = {-alpha_p4 / policy.set_index('coordinate').loc['k', 'port_dRe_s']:.4f} vs actual {bis_res['k'] - P4[1]:.4f}",
            "newton": f"{len(k_newton) - 1} steps to k = {k_newton[-1]['k']:.10f} (|Re s| = {abs(k_newton[-1]['re_s']):.1e})",
            "frozen_reference": "F7A grid bracket",
            "tds": "none",
            "independent_tool": "no (custom GFL)",
        },
        {
            "intervention": "damped synchronous condenser (G1)",
            "initial_alpha": alpha_p4,
            "final_state": "; ".join(
                f"{r.condenser:.4f}: {r.alpha_flag_frozen:+.4f}"
                for r in cond.itertuples()
                if r.condenser == r.condenser
            ),
            "port_local_derivative": "N/A (device outside the port family)",
            "one_step_prediction": "N/A",
            "newton": "N/A",
            "frozen_reference": "threshold 2.48 % = 106.0 MVA",
            "tds": "frozen G2: 2 % unstable (+0.0264), 3 % stable",
            "independent_tool": "no",
        },
        {
            "intervention": "documented TGOV1N governors (FC03)",
            "initial_alpha": alpha_p4,
            "final_state": f"alpha = {gov.alpha_flag_governed:+.4f}; H governed = {gov.H_governed}",
            "port_local_derivative": "N/A (model change)",
            "one_step_prediction": "N/A",
            "newton": "N/A",
            "frozen_reference": "FC03",
            "tds": "none",
            "independent_tool": "governor damping direction only (ANDES, FC14)",
        },
        {
            "intervention": "single-line reinforcement (12 frozen lines)",
            "initial_alpha": alpha_p4,
            "final_state": f"max stabilization at 2x: {-(lines_df['actual_delta_alpha_2.0']).min():.4f} s^-1 ({lines_df.loc[lines_df['actual_delta_alpha_2.0'].idxmin(), 'line']}); any line restores stability at <= 2x: {bool((lines_df['status_2.0'] == 'STABLE').any())}",
            "port_local_derivative": "per line (PCV04_lines.csv)",
            "one_step_prediction": f"Spearman port vs actual (1.5x) = {task_e.iloc[0]['spearman']:.3f}",
            "newton": "N/A",
            "frozen_reference": "F2c holdout (boundary derivatives, GFL)",
            "tds": "none",
            "independent_tool": "equation-equivalent R3 line sensitivities reproduced in ANDES (0a4e18c2), different converter model",
        },
    ]
    remediation = pd.DataFrame(rem)
    remediation.to_csv(OUT / "PCV04_remediation.csv", index=False)
    write_json(
        OUT / "PCV04_summary.json",
        {
            "alpha_p4": alpha_p4,
            "lambda_p4": [lam_p4.real, lam_p4.imag],
            "s_star_c2": [s_star.real, s_star.imag],
            "newton_g": g_newton,
            "newton_k": k_newton,
            "bisect": bis_res,
            "gscr0": gscr0,
        },
    )
    pd.set_option("display.width", 250)
    print(policy.to_string(index=False))
    print(lines_df.to_string(index=False))
    print(task_e.to_string(index=False))
    print(remediation.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
