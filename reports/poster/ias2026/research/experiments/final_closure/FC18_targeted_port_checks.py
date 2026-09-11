"""FC18: targeted port-structure checks on the frozen IEEE-39 (user's final transfer).

Only the frozen model is used (core 30/33/35/37, matched dispatch, no GFM, no DC
link, no limits, no governor). Objects:

    descriptor pencil  P_S(s) = [[sI - f_x, -f_z], [g_x, g_z]]  on the common
                       realization C2 (BC02: every candidate carries both devices;
                       the Jacobian is affine in the binary indicators)
    port               T_S(s) = g_z + g_x (sI - f_x)^-1 f_z   (network variables kept)
    small port matrix  M(s) = D(s) K(s),  D = blkdiag(dY_i),  K = E^T T_0^-1 E
                       det T_S / det T_0 = det(I + M_SS)
    interaction        I + M = (I + D K_d)(I + Q),  Q = (I + D K_d)^-1 D K_o,
                       K_d / K_o the diagonal / off-diagonal 2x2 blocks of K
                       (network closure between candidate buses lives in K_o)

Checks
    1  descriptor affinity per block; ratio identity for the 16 subsets; principal-
       minor (Moebius) expansion by order at the critical inter-area boundary
    2  rank / singular spectrum of Q along the frozen F7B path (t = 0.8515625) whose
       minimum witness goes 4 -> 3 -> 2 -> 3 -> 4
    3  boundary gradient ds*/dtheta = -(p* F_theta q)/(p* F_s q), F = I + M_SS,
       against direct full-DAE finite differences; Q/V retuning prediction at P4
    4  walk expansion of T_0^-1 (bus-diagonal splitting); NON-EXECUTABLE if the
       Neumann condition rho(D_g^-1 N) < 1 fails
(5, the second-order spot check, is FC07 v2; 6 is a rule, not a computation.)
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _fc import WORKERS, FCExperiment, out_dir, write_json
from BC02_binary_representation import port_t  # noqa: E402

from _f7_common import (  # noqa: E402
    LEAK,
    Theta,
    native_te,
    solve_subset,
    te_geometric_mean,
)
from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    all_vertices,
    build_common_realization,
    reduced,
)
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402

OUT = out_dir("FC18_targeted_port_checks")
CORE = (30, 33, 35, 37)
PATH_T = 0.8515625  # frozen F7B grid line (k = 1, h = 1)
# (label, subset whose label changes, g_lo, g_hi) from the frozen F7B grid
EVENTS = [
    ("E1 kappa 4->3", (30, 33, 35), 0.01, 0.0103125),
    ("E2", (30, 33, 37), 0.01875, 0.0190625),
    ("E3", (33, 35, 37), 0.060625, 0.06125),
    ("E4 kappa 3->2", (30, 33), 0.103125, 0.10375),
    ("E5", (33, 35, 37), 0.186875, 0.1878125),
    ("E6 kappa 2->3", (30, 33), 0.2484375, 0.25),
    ("E7", (30, 33, 35), 0.2765625, 0.278125),
    ("E8 kappa 3->4", (30, 33, 37), 0.2921875, 0.29375),
    ("E9 -> EMPTY", (30, 33, 35, 37), 0.296875, 0.2984375),
]
P4 = (0.03625, 1.425, 1.5, 1.0)
FLAG = (
    "FLAG P4-line",
    CORE,
)  # boundary along g at k = 1.425, t = 1.5 (FC13 g* = 0.2077)
STEPS = {"g": 1e-4, "k": 1e-4, "t": 1e-4, "h": 1e-4}
NAMES = ("g", "k", "t", "h")


# ------------------------------------------------------------------ models --
def kwargs_of(th):
    g, k, t, h = th
    kw = {
        "converter": ConverterParameters(
            voltage_control=True, voltage_gain=g, voltage_leak=LEAK
        ),
        "machine_scaling": {"ka": k, "ta": t},
    }
    if h != 1.0:
        mean = te_geometric_mean()
        kw["machine_bus_scaling"] = {
            b: {"ta": (mean / te) ** (1.0 - h)} for b, te in native_te().items()
        }
    return kw


def direct_perp(members, th):
    case = solve_subset(tuple(members), Theta(*th))
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(case.system.A, r_x, frequency_partner(case.dae).w)
    return np.linalg.eigvals(tr.a_perp)


def critical(ev):
    """The rightmost oscillatory transverse eigenvalue (upper half plane)."""

    osc = ev[ev.imag > 0.5]
    return osc[np.argmax(osc.real)]


def locate(members, th_fn, lo, hi, iters=40):
    a_lo = critical(direct_perp(members, th_fn(lo))).real
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        a = critical(direct_perp(members, th_fn(mid))).real
        if np.sign(a) == np.sign(a_lo):
            lo, a_lo = mid, a
        else:
            hi = mid
    g = 0.5 * (lo + hi)
    return g, critical(direct_perp(members, th_fn(g)))


class Realization:
    def __init__(self, members, th=None):
        self.members = tuple(members)
        kw = {} if th is None else kwargs_of(th)
        self.cr = build_common_realization(self.members, **kw)
        self.jac = {}
        for delta, s in all_vertices(self.members):
            j, _, _ = self.cr.jacobian(delta, "C2")
            self.jac[tuple(sorted(s))] = j
        net = self.cr.case.dae.network
        self.ch = {
            b: [2 * net.position(b), 2 * net.position(b) + 1] for b in self.members
        }
        self.idx = [c for b in self.members for c in self.ch[b]]

    def ports(self, s):
        return {t: port_t(j, s) for t, j in self.jac.items()}

    def m_matrix(self, s):
        p = self.ports(s)
        t0 = p[()]
        kfull = np.linalg.inv(t0)
        k = kfull[np.ix_(self.idx, self.idx)]
        n = 2 * len(self.members)
        d = np.zeros((n, n), complex)
        for i, b in enumerate(self.members):
            c = self.ch[b]
            d[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (p[(b,)] - t0)[np.ix_(c, c)]
        return d, k, p, t0

    def cols(self, subset):
        return [2 * self.members.index(b) + c for b in subset for c in (0, 1)]

    def vertex_eig(self, subset, target):
        ev = np.linalg.eigvals(reduced(self.jac[tuple(sorted(subset))]))
        return ev[np.argmin(np.abs(ev - target))]


def subsets(items):
    for r in range(len(items) + 1):
        yield from combinations(items, r)


def moebius(values, items):
    return {
        t: sum((-1) ** (len(t) - len(q)) * values[q] for q in subsets(t))
        for t in subsets(items)
    }


def split_interaction(d, k):
    n = k.shape[0]
    kd = np.zeros_like(k)
    for i in range(0, n, 2):
        kd[i : i + 2, i : i + 2] = k[i : i + 2, i : i + 2]
    ko = k - kd
    loc = np.eye(n) + d @ kd
    q = np.linalg.solve(loc, d @ ko)
    return kd, ko, loc, q


# ----------------------------------------------------------------- check 1 --
def check1(r: Realization, subset, s):
    """Affinity per block, ratio identity (16 subsets), expansions at s."""

    items = r.members
    blocks = {}
    for name in ("fx", "fz", "gx", "gz"):
        vals = {t: getattr(r.jac[t], name) for t in r.jac}
        mu = moebius(vals, items)
        base = np.linalg.norm(vals[()])
        by = {}
        for t, c in mu.items():
            by[len(t)] = max(by.get(len(t), 0.0), float(np.linalg.norm(c) / base))
        support = [int((np.abs(mu[(b,)]) > 1e-12 * base).sum()) for b in items]
        blocks[name] = {"moebius_rel_by_order": by, "singleton_increment_nnz": support}
    d, k, p, t0 = r.m_matrix(s)
    m = d @ k
    s0, l0 = np.linalg.slogdet(t0)
    ident = 0.0
    ratio, det_ipm = {}, {}
    for t in subsets(items):
        st, lt = np.linalg.slogdet(p[t])
        ratio[t] = (st / s0) * np.exp(lt - l0)
        c = r.cols(t)
        det_ipm[t] = np.linalg.det(np.eye(len(c)) + m[np.ix_(c, c)]) if c else 1.0 + 0j
        ident = max(ident, abs(ratio[t] - det_ipm[t]) / max(1.0, abs(det_ipm[t])))
    sub_items = tuple(sorted(subset))
    mu = moebius({t: det_ipm[t] for t in subsets(sub_items)}, sub_items)
    by_order = {
        kk: complex(sum(v for t, v in mu.items() if len(t) == kk))
        for kk in range(len(sub_items) + 1)
    }
    trunc = {kk: abs(sum(by_order[j] for j in range(kk + 1))) for kk in by_order}
    kd, ko, loc, q = split_interaction(d, k)
    c = r.cols(sub_items)
    local = {
        t: np.linalg.det(
            np.eye(len(r.cols(t))) + (d @ kd)[np.ix_(r.cols(t), r.cols(t))]
        )
        if t
        else 1.0 + 0j
        for t in subsets(sub_items)
    }
    inter = {
        t: np.linalg.det(np.eye(len(r.cols(t))) + q[np.ix_(r.cols(t), r.cols(t))])
        if t
        else 1.0 + 0j
        for t in subsets(sub_items)
    }
    mu_q = moebius(inter, sub_items)
    q_order = {
        kk: complex(sum(v for t, v in mu_q.items() if len(t) == kk))
        for kk in range(len(sub_items) + 1)
    }
    q_trunc = {kk: abs(sum(q_order[j] for j in range(kk + 1))) for kk in q_order}
    fact_err = abs(det_ipm[sub_items] - local[sub_items] * inter[sub_items])
    return {
        "blocks": blocks,
        "ratio_identity_max_rel": ident,
        "det_I_plus_M_S": abs(det_ipm[sub_items]),
        "r_order_sums": {kk: [v.real, v.imag, abs(v)] for kk, v in by_order.items()},
        "r_truncation_abs": trunc,
        "local_factor_abs": abs(local[sub_items]),
        "interaction_factor_abs": abs(inter[sub_items]),
        "interaction_order_sums": {
            kk: [v.real, v.imag, abs(v)] for kk, v in q_order.items()
        },
        "interaction_truncation_abs": q_trunc,
        "factorization_error": fact_err,
        "local_only_model_det": abs(local[sub_items]),
    }, (d, k, q)


# ----------------------------------------------------------------- check 2 --
def check2(d, k, q, r: Realization, s):
    sv_q = np.linalg.svd(q, compute_uv=False)
    kd, ko, _, _ = split_interaction(d, k)
    sv_ko = np.linalg.svd(ko, compute_uv=False)
    sv_m = np.linalg.svd(d @ k, compute_uv=False)
    inter = {
        t: np.linalg.det(np.eye(len(r.cols(t))) + q[np.ix_(r.cols(t), r.cols(t))])
        if t
        else 1.0 + 0j
        for t in subsets(r.members)
    }
    mu = moebius(inter, r.members)
    order = {
        kk: max(abs(v) for t, v in mu.items() if len(t) == kk) for kk in range(2, 5)
    }
    return {
        "sv_Q": sv_q.tolist(),
        "sv_Ko": sv_ko.tolist(),
        "sv_M": sv_m.tolist(),
        "eff_rank_Q_1e-1": int((sv_q / sv_q[0] > 1e-1).sum()),
        "eff_rank_Q_1e-2": int((sv_q / sv_q[0] > 1e-2).sum()),
        "eff_rank_Q_1e-3": int((sv_q / sv_q[0] > 1e-3).sum()),
        "max_abs_interaction_coeff_by_order": order,
        "order4_over_prod_top4_sv": float(order[4] / np.prod(sv_q[:4])),
    }


# ----------------------------------------------------------------- check 3 --
def crit_pair(m_ss):
    ev, vr = np.linalg.eig(m_ss)
    i = int(np.argmin(np.abs(ev + 1.0)))
    evl, vl = np.linalg.eig(m_ss.conj().T)
    j = int(np.argmin(np.abs(evl - np.conj(ev[i]))))
    return ev[i], vr[:, i], vl[:, j]


def port_gradient(subset, th, s_star, center: Realization):
    c = center.cols(subset)

    def mss(r, s):
        d, k, _, _ = r.m_matrix(s)
        return (d @ k)[np.ix_(c, c)]

    m0 = mss(center, s_star)
    lam, q, p = crit_pair(m0)
    hs = 1e-5
    f_s = (mss(center, s_star + hs) - mss(center, s_star - hs)) / (2 * hs)
    den = p.conj() @ f_s @ q
    grad = {}
    for i, name in enumerate(NAMES):
        dp, dm = list(th), list(th)
        dp[i] += STEPS[name]
        dm[i] -= STEPS[name]
        f_th = (
            mss(Realization(CORE, tuple(dp)), s_star)
            - mss(Realization(CORE, tuple(dm)), s_star)
        ) / (2 * STEPS[name])
        grad[name] = complex(-(p.conj() @ f_th @ q) / den)
    return grad, complex(lam)


def direct_gradient(subset, th, s_star):
    grad = {}
    for i, name in enumerate(NAMES):
        dp, dm = list(th), list(th)
        dp[i] += STEPS[name]
        dm[i] -= STEPS[name]
        ep = direct_perp(subset, dp)
        em = direct_perp(subset, dm)
        lp = ep[np.argmin(np.abs(ep - s_star))]
        lm = em[np.argmin(np.abs(em - s_star))]
        grad[name] = complex((lp - lm) / (2 * STEPS[name]))
    return grad


SCALE = np.array([1.0, 1.8, 2.5, 2.0])  # policy-box ranges (g, k, t, h)


def compare(gp, gd):
    vp = np.array([gp[n].real for n in NAMES]) * SCALE
    vd = np.array([gd[n].real for n in NAMES]) * SCALE
    cos = float(vp @ vd / (np.linalg.norm(vp) * np.linalg.norm(vd)))
    ip = np.array([gp[n].imag for n in NAMES]) * SCALE
    idd = np.array([gd[n].imag for n in NAMES]) * SCALE
    return {
        "normal_angle_deg": float(np.degrees(np.arccos(np.clip(cos, -1, 1)))),
        "normal_magnitude_rel_err": float(
            abs(np.linalg.norm(vp) - np.linalg.norm(vd)) / np.linalg.norm(vd)
        ),
        "dRe_dg_port": gp["g"].real,
        "dRe_dg_direct": gd["g"].real,
        "dIm_rel_err": float(
            np.linalg.norm(ip - idd) / max(np.linalg.norm(idd), 1e-300)
        ),
    }


# ------------------------------------------------------------------ events --
def event_task(ev):
    pin_blas_threads()
    label, subset, lo, hi, th_fn_kind = ev
    if th_fn_kind == "path":

        def th_fn(g):
            return (g, 1.0, PATH_T, 1.0)
    else:

        def th_fn(g):
            return (g, 1.425, 1.5, 1.0)

    g_star, lam_direct = locate(subset, th_fn, lo, hi)
    th = th_fn(g_star)
    r = Realization(CORE, th)
    s_star = r.vertex_eig(subset, lam_direct)
    c1, (d, k, q) = check1(r, subset, s_star)
    c2 = check2(d, k, q, r, s_star)
    gp, lam_m = port_gradient(subset, th, s_star, r)
    gd = direct_gradient(subset, th, lam_direct)
    return {
        "event": label,
        "subset": "+".join(map(str, subset)),
        "g_star": g_star,
        "theta": th,
        "lambda_direct": [lam_direct.real, lam_direct.imag],
        "s_star_realization": [s_star.real, s_star.imag],
        "freq_hz": float(abs(s_star.imag) / (2 * np.pi)),
        "critical_M_eig_plus_1": abs(lam_m + 1.0),
        "check1": c1,
        "check2": c2,
        "grad_port": {n: [v.real, v.imag] for n, v in gp.items()},
        "grad_direct": {n: [v.real, v.imag] for n, v in gd.items()},
        "check3": compare(gp, gd),
    }


# ------------------------------------------------------ retuning at P4 --
def retune_task(_):
    """Predict the minimum Q/V gain change that makes the P4 flagship stable."""

    pin_blas_threads()
    th = list(P4)
    iters = []
    for it in range(6):
        lam = critical(direct_perp(CORE, th))
        r = Realization(CORE, tuple(th))
        s = r.vertex_eig(CORE, lam)
        c = r.cols(CORE)

        def mss(rr, ss, c=c):
            d, k, _, _ = rr.m_matrix(ss)
            return (d @ k)[np.ix_(c, c)]

        _, qv, pv = crit_pair(mss(r, s))
        hs = 1e-5
        f_s = (mss(r, s + hs) - mss(r, s - hs)) / (2 * hs)
        dp, dm = list(th), list(th)
        dp[0] += STEPS["g"]
        dm[0] -= STEPS["g"]
        f_g = (
            mss(Realization(CORE, tuple(dp)), s) - mss(Realization(CORE, tuple(dm)), s)
        ) / (2 * STEPS["g"])
        dsdg = complex(-(pv.conj() @ f_g @ qv) / (pv.conj() @ f_s @ qv))
        step = -s.real / dsdg.real
        iters.append(
            {
                "iter": it,
                "g": th[0],
                "re_s": s.real,
                "im_s": s.imag,
                "dRe_dg_port": dsdg.real,
                "step": step,
            }
        )
        if abs(s.real) < 1e-7:
            break
        th[0] += step
    return iters


# ------------------------------------------------------------ check 4 --
def walk_check(r: Realization, subset, s, label):
    d, k, p, t0 = r.m_matrix(s)
    n = t0.shape[0]
    dg = np.zeros_like(t0)
    for i in range(0, n, 2):
        dg[i : i + 2, i : i + 2] = t0[i : i + 2, i : i + 2]
    nmat = dg - t0
    wmat = np.linalg.solve(dg, nmat)
    rho = float(np.abs(np.linalg.eigvals(wmat)).max())
    out = {"case": label, "s": [s.real, s.imag], "rho_walk": rho}
    if rho >= 1.0:
        out["status"] = "NON-EXECUTABLE (Neumann condition fails: rho(D_g^-1 N) >= 1)"
        # the diagnostic that is still meaningful: how much of the critical coupling
        # is network closure (off-diagonal K) rather than local (diagonal K)
        c = r.cols(subset)
        m = d @ k
        lam, q, pv = crit_pair(m[np.ix_(c, c)])
        kd, ko, _, _ = split_interaction(d, k)
        dsub = d[np.ix_(c, c)]
        tot = pv.conj() @ m[np.ix_(c, c)] @ q
        net = pv.conj() @ (dsub @ ko[np.ix_(c, c)]) @ q
        out["critical_coupling_network_share"] = float(abs(net) / abs(tot))
        return out
    # convergent: walk-length series and path contributions
    dinv = np.linalg.inv(dg)
    kfull = np.linalg.inv(t0)
    acc = dinv.copy()
    term = dinv.copy()
    conv = []
    for _length in range(1, 61):
        term = wmat @ term
        acc = acc + term
        conv.append(
            float(
                np.linalg.norm((acc - kfull)[np.ix_(r.idx, r.idx)])
                / np.linalg.norm(kfull[np.ix_(r.idx, r.idx)])
            )
        )
    out["status"] = "CONVERGENT"
    out["series_rel_error_by_length"] = conv[::5]
    return out


def control_task(_):
    pin_blas_threads()
    members = (31, 32, 33, 38)
    case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    ev = np.linalg.eigvals(
        transverse_operator(case.system.A, r_x, frequency_partner(case.dae).w).a_perp
    )
    band = ev[(ev.imag > 2 * np.pi * 0.3) & (ev.imag < 2 * np.pi * 1.5)]
    lam = band[np.argmax(band.real)]
    r = Realization(members, None)
    s = r.vertex_eig(members, lam)
    return walk_check(
        r, members, s, "matched control 31+32+33+38 (E12 census policy)"
    ), [lam.real, lam.imag]


def flag_walk_task(ev_result):
    pin_blas_threads()
    th = tuple(ev_result["theta"])
    r = Realization(CORE, th)
    s = complex(*ev_result["s_star_realization"])
    return walk_check(r, CORE, s, "flagship boundary (P4 line)")


def main(argv) -> int:
    exp = FCExperiment(
        name="FC18_targeted_port_checks",
        question=(
            "Do descriptor closure, interaction rank, port gradients and network "
            "walks explain the portfolio structure?"
        ),
        config={"path_t": PATH_T, "events": [e[0] for e in EVENTS], "steps": STEPS},
        workers=WORKERS,
    )
    started = time.time()
    tasks = [(lab, sub, lo, hi, "path") for lab, sub, lo, hi in EVENTS]
    tasks.append((FLAG[0], FLAG[1], 0.20, 0.22, "p4line"))
    with Pool(min(WORKERS, len(tasks) + 2)) as pool:
        ev_async = pool.map_async(event_task, tasks, chunksize=1)
        rt_async = pool.map_async(retune_task, [0])
        ct_async = pool.map_async(control_task, [0])
        events = ev_async.get()
        retune = rt_async.get()[0]
        control, control_lam = ct_async.get()[0]
    flag = next(e for e in events if e["event"] == FLAG[0])
    flag_walk = flag_walk_task(flag)
    summary = {
        "events": events,
        "retune_P4": retune,
        "retune_direct_boundary_g": 0.20768140450381395,
        "walk": [flag_walk, control],
        "control_lambda": control_lam,
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC18_summary.json", summary)
    rows = []
    for e in events:
        c1, c2, c3 = e["check1"], e["check2"], e["check3"]
        rows.append(
            {
                "event": e["event"],
                "subset": e["subset"],
                "g_star": e["g_star"],
                "freq_hz": e["freq_hz"],
                "ratio_identity": c1["ratio_identity_max_rel"],
                "det_I_plus_M": c1["det_I_plus_M_S"],
                "local_factor": c1["local_factor_abs"],
                "interaction_factor": c1["interaction_factor_abs"],
                "trunc_interaction": json.dumps(
                    {
                        k: round(v, 4)
                        for k, v in c1["interaction_truncation_abs"].items()
                    }
                ),
                "sv_Q": json.dumps([round(x, 4) for x in c2["sv_Q"]]),
                "eff_rank_Q_1e-2": c2["eff_rank_Q_1e-2"],
                "order4_coeff": c2["max_abs_interaction_coeff_by_order"][4],
                "normal_angle_deg": c3["normal_angle_deg"],
                "normal_mag_err": c3["normal_magnitude_rel_err"],
                "dIm_rel_err": c3["dIm_rel_err"],
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "FC18_events.csv", index=False)
    exp.finish("COMPUTED", elapsed_s=summary["elapsed_s"])
    print(pd.DataFrame(rows).to_string(index=False))
    print(
        json.dumps({"retune": retune, "walk": summary["walk"]}, indent=1, default=str)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
