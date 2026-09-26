# ruff: noqa: E501
"""E11 spectral graph baselines and E12 controller-weighted modal mixing (prereg E11, E12)."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

PHASE12 = "E12"
PHASE12Q3 = "E12q3"
POLICIES = list(C.DISCOVERY) + [f"H{i:02d}" for i in range(1, 25)]
CORE16 = C.subsets(C.V4)


def graph_basis():
    L, order = C.laplacian_B()
    w, U = np.linalg.eigh(L)
    return w, U, order


def angle_shape(phi_rect, V0):
    """Bus-voltage angle perturbation from a rectangular mode shape (complex 78-vector)."""

    dv = phi_rect[0::2] + 1j * phi_rect[1::2]  # complex bus-voltage perturbation per bus
    return np.imag(np.conj(V0) * dv) / np.abs(V0) ** 2, np.real(np.conj(V0) * dv) / np.abs(V0)


def graph_energy(vec, U):
    c = U.T @ vec
    e = np.abs(c) ** 2
    return e / e.sum() if e.sum() > 0 else e


def support80(energy):
    idx = np.argsort(-energy)
    cum = np.cumsum(energy[idx])
    k = int(np.searchsorted(cum, 0.8) + 1)
    return frozenset(int(i) for i in idx[:k])


# ------------------------------------------------------------------ E12 tasks --
def tasks12():
    out = [{"phase": PHASE12, "pid": p, "S": list(s)} for p in POLICIES for s in CORE16]
    return out


def mixing(T_s, U):
    K = np.kron(U, np.eye(2))
    That = K.T @ T_s @ K
    B = np.zeros_like(That)
    for k in range(U.shape[1]):
        B[2 * k: 2 * k + 2, 2 * k: 2 * k + 2] = That[2 * k: 2 * k + 2, 2 * k: 2 * k + 2]
    R = That - B
    return float(np.linalg.norm(R) / np.linalg.norm(That))


def run_task12(task):
    theta = C.policy(task["pid"])
    try:
        case = C.solve(tuple(task["S"]), theta)
    except (C.InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as e:
        return {"status": "INFEASIBLE", "error": str(e)[:100]}
    a, d, jac = C.physical_matrices(case)
    ev = C.transverse_eval(case, modes=False)
    w, U, order = graph_basis()
    s = 1j * 2 * np.pi * ev["lam_hz"] if ev["lam_hz"] > 1e-3 else 1j * 2 * np.pi * 0.7
    n = jac.fx.shape[0]
    T = jac.gz + jac.gx @ np.linalg.solve(s * np.eye(n) - jac.fx, jac.fz)
    mu = mixing(T, U)
    mu_net = mixing(jac.gz, U)  # algebraic part only (network + loads + device algebraic terms)
    # critical-mode graph-Fourier energy (angle shape)
    vals, vecs = np.linalg.eig(a)
    lam_c = complex(ev["lam_re"], 2 * np.pi * ev["lam_hz"])
    j = int(np.argmin(np.abs(vals - lam_c)))
    phi = -np.linalg.solve(jac.gz, jac.gx @ vecs[:, j])
    V0 = case.equilibrium.z[0::2] + 1j * case.equilibrium.z[1::2]
    ang, mag = angle_shape(phi, V0)
    e_ang = graph_energy(ang, U)
    return {"status": ev["status"], "alpha": ev["alpha"], "lam_hz": ev["lam_hz"], "mu_mix": mu, "mu_mix_gz": mu_net,
            "energy_angle": e_ang.round(6).tolist(), "low3_share": float(e_ang[1:4].sum()), "top3_modes": np.argsort(-e_ang)[:3].tolist()}


# -------------------------------------------------------------- E12 Q3 tasks --
def fc18_events():
    data = json.loads((C.RESEARCH / "outputs/ias2026/final_math_nonlinear_validation_20260910T231539/FC18_targeted_port_checks/FC18_summary.json").read_text())
    ev = [e for e in data["events"] if e["event"].startswith("E")]
    return sorted(({"event": e["event"], "g": e["theta"][0], "theta": e["theta"]} for e in ev), key=lambda x: x["g"])


def tasks12q3():
    evs = fc18_events()
    gs = [e["g"] for e in evs]
    ctrl = [0.5 * gs[0]] + [0.5 * (gs[k] + gs[k + 1]) for k in range(len(gs) - 1)]
    out = []
    th = evs[0]["theta"]
    for e in evs:
        for side, f in (("minus", 0.99), ("plus", 1.01)):
            out.append({"phase": PHASE12Q3, "kind": "event", "name": e["event"], "side": side, "theta": [e["g"] * f] + list(th[1:])})
    for k, gc in enumerate(ctrl):
        for side, f in (("minus", 0.99), ("plus", 1.01)):
            out.append({"phase": PHASE12Q3, "kind": "control", "name": f"C{k}", "side": side, "theta": [gc * f] + list(th[1:])})
    return out


def run_task12q3(task):
    theta = tuple(task["theta"])
    w, U, order = graph_basis()
    stat, alpha, sup = {}, {}, {}
    for s in CORE16:
        case = C.solve(s, theta)
        ev = C.transverse_eval(case, modes=False)
        stat[C.label(s)], alpha[C.label(s)] = ev["status"], ev["alpha"]
    unsafe = [s for s in CORE16 if stat[C.label(s)] == "UNSTABLE"]
    minimal = [s for s in unsafe if not any(set(r) < set(s) for r in unsafe)]
    H = "|".join(sorted(C.label(e) for e in minimal)) or "EMPTY"
    out = {"H": H}
    if minimal:
        wit = max(minimal, key=lambda s: alpha[C.label(s)])
        case = C.solve(wit, theta)
        a, _, jac = C.physical_matrices(case)
        ev = C.transverse_eval(case, modes=False)
        vals, vecs = np.linalg.eig(a)
        j = int(np.argmin(np.abs(vals - complex(ev["lam_re"], 2 * np.pi * ev["lam_hz"]))))
        phi = -np.linalg.solve(jac.gz, jac.gx @ vecs[:, j])
        V0 = case.equilibrium.z[0::2] + 1j * case.equilibrium.z[1::2]
        e = graph_energy(angle_shape(phi, V0)[0], U)
        sup = sorted(support80(e))
        out.update(witness=C.label(wit), support=sup, lam_hz=ev["lam_hz"])
    return out


# ---------------------------------------------------------------- aggregates --
def aggregate12():
    rows = []
    for r in I.Store(PHASE12).all():
        if not r.get("ok") or r.get("status") == "INFEASIBLE":
            continue
        t = r["task"]
        rows.append({"pid": t["pid"], "S": C.label(t["S"]), **{k: r[k] for k in ("status", "alpha", "lam_hz", "mu_mix", "mu_mix_gz", "low3_share")},
                     "top3_modes": "-".join(map(str, r["top3_modes"]))})
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E12_mixing.csv", index=False)
    h4 = df[df.S == "30+33+35+37"].set_index("pid")
    q1 = float((h4.mu_mix.max() - h4.mu_mix.min()) / h4.mu_mix.median())
    out = {"Q1_rel_range_mu_H4": q1, "Q1_material": bool(q1 >= 0.20), "mu_mix_H4_min": float(h4.mu_mix.min()),
           "mu_mix_H4_max": float(h4.mu_mix.max()), "mu_mix_gz_H4_rel_range": float((h4.mu_mix_gz.max() - h4.mu_mix_gz.min()) / h4.mu_mix_gz.median())}
    pol = pd.read_csv(I.RESULTS / "CDW_E1_policy_summary.csv") if (I.RESULTS / "CDW_E1_policy_summary.csv").exists() else None
    if pol is not None:
        hp = pol[(pol.split == "holdout") & pol.base_stable].set_index("pid")
        x = h4.mu_mix.reindex(hp.index).to_numpy(float)
        y = hp.n_rev_interventions.to_numpy(float)
        rho = AN.spearman(x, y)
        rng = np.random.default_rng(20260924)
        ok = np.isfinite(x) & np.isfinite(y)
        perm = [AN.spearman(x[ok], rng.permutation(y[ok])) for _ in range(10000)]
        p = float((np.abs(perm) >= abs(rho)).mean()) if np.isfinite(rho) else np.nan
        out.update(Q2_spearman=rho, Q2_perm_p=p, Q2_correlates=bool(np.isfinite(rho) and abs(rho) >= 0.5 and p <= 0.01))
    # Q4: top-3 dynamically critical graph modes of H4 at D01 vs 3 lowest nonzero L_B modes
    t = h4.loc["D01", "top3_modes"] if "D01" in h4.index else ""
    top3 = {int(v) for v in t.split("-")} if t else set()
    low3 = {1, 2, 3}
    out["Q4_top3_dynamic"] = sorted(top3)
    out["Q4_jaccard_vs_low3"] = len(top3 & low3) / len(top3 | low3) if top3 else np.nan
    # Q3
    q3 = {}
    for r in I.Store(PHASE12Q3).all():
        if r.get("ok"):
            t = r["task"]
            q3.setdefault((t["kind"], t["name"]), {})[t["side"]] = r
    ev_changes = ctrl_changes = n_ev = n_ctrl = 0
    q3rows = []
    for (kind, name), sides in q3.items():
        if len(sides) < 2:
            continue
        a, b = sides["minus"], sides["plus"]
        wchg = a.get("H") != b.get("H")
        schg = a.get("support") != b.get("support")
        q3rows.append({"kind": kind, "name": name, "H_minus": a.get("H"), "H_plus": b.get("H"), "witness_changed": wchg,
                       "support_minus": a.get("support"), "support_plus": b.get("support"), "support_changed": schg})
        if kind == "event":
            n_ev += 1
            ev_changes += schg
        else:
            n_ctrl += 1
            ctrl_changes += schg
    pd.DataFrame(q3rows).to_csv(I.RESULTS / "CDW_E12_q3_support.csv", index=False)
    out.update(Q3_events_with_support_change=int(ev_changes), Q3_n_events=int(n_ev), Q3_controls_with_support_change=int(ctrl_changes),
               Q3_n_controls=int(n_ctrl), Q3_reproducible=bool(ev_changes >= 7 and n_ctrl and ctrl_changes / n_ctrl <= 0.20))
    I.atomic_write_json(I.RESULTS / "CDW_E12_gates.json", out)
    return out


def aggregate11():
    """Static graph scores: bases, and predictive tests against E1/E3/E4/E10 data (no new solves)."""

    w, U, order = graph_basis()
    idx = {b: i for i, b in enumerate(order)}
    L, _ = C.laplacian_B()
    Lp = np.linalg.pinv(L)
    u2 = U[:, 1]
    node = []
    for b in C.V9:
        e = np.zeros(len(order))
        e[idx[b]], e[idx[39]] = 1, -1
        node.append({"bus": b, "fiedler_entry": float(u2[idx[b]]), "abs_fiedler": abs(float(u2[idx[b]])), "reff_to_39": float(e @ Lp @ e)})
    nd = pd.DataFrame(node).set_index("bus")
    # Kron reduction of B_bus onto the 10 generator buses (basis B)
    Y = C.build_ybus()
    B = -Y.imag
    gen = [idx[b] for b in C.SG_BUSES]
    oth = [i for i in range(len(order)) if i not in gen]
    Bk = B[np.ix_(gen, gen)] - B[np.ix_(gen, oth)] @ np.linalg.solve(B[np.ix_(oth, oth)], B[np.ix_(oth, gen)])
    Lk = -np.triu(Bk, 1) - np.triu(Bk, 1).T
    Lk = Lk - np.diag(Lk.sum(1))
    Lk = -Lk if np.trace(Lk) < 0 else Lk
    wk = np.linalg.eigvalsh(Lk)
    out = {"LB_eigs_low5": w[:5].round(6).tolist(), "kron_eigs": wk.round(6).tolist()}
    tests = []
    e10 = pd.read_csv(I.RESULTS / "CDW_E10_shapley.csv") if (I.RESULTS / "CDW_E10_shapley.csv").exists() else None
    e1 = pd.read_csv(I.RESULTS / "CDW_E1_summary.csv") if (I.RESULTS / "CDW_E1_summary.csv").exists() else None
    e3 = pd.read_parquet(I.RESULTS / "CDW_E3_node_sensitivity.parquet") if (I.RESULTS / "CDW_E3_node_sensitivity.parquet").exists() else None
    for score in ("fiedler_entry", "abs_fiedler", "reff_to_39"):
        s = nd[score]
        for pid in [f"H{i:02d}" for i in range(1, 25)]:
            if e1 is not None:
                g = e1[e1.pid == pid].set_index("i")
                tests.append({"score": score, "target": "E1_frac_destab", "pid": pid, "rho": AN.spearman(s.reindex(g.index), g.frac_destab)})
            if e10 is not None:
                g = e10[e10.pid == pid].set_index("i")
                if len(g):
                    tests.append({"score": score, "target": "E10_shapley", "pid": pid, "rho": AN.spearman(s.reindex(g.index), g.shapley)})
            if e3 is not None:
                g = e3[(e3.pid == pid) & (e3.target == "V9") & (e3.kind == "g") & (e3.semantics == "SPR")].set_index("idx")
                if len(g):
                    tests.append({"score": score, "target": "E3_NG_total", "pid": pid, "rho": AN.spearman(s.reindex(g.index), g.d_total)})
    td = pd.DataFrame(tests)
    td.to_csv(I.RESULTS / "CDW_E11_spectral.csv", index=False)
    summ = td.groupby(["score", "target"]).rho.apply(lambda x: float(np.nanmedian(np.abs(x)))).rename("median_abs_rho").reset_index()
    out["node_tests"] = summ.to_dict("records")
    out["static_graph_predicts_any"] = bool((summ.median_abs_rho >= 0.6).any())
    bl = pd.read_csv(I.RESULTS / "CDW_E4_baselines.csv") if (I.RESULTS / "CDW_E4_baselines.csv").exists() else None
    if bl is not None:
        h = bl[bl.split.str.startswith("holdout") & (bl.target == "H4")]
        out["edge_S5_median_rho"] = float(h[h.predictor == "S5_reff"].spearman.median())
        out["edge_S6_median_rho"] = float(h[h.predictor == "S6_fiedler"].spearman.median())
    nd.to_csv(I.RESULTS / "CDW_E11_node_graph_scores.csv")
    I.atomic_write_json(I.RESULTS / "CDW_E11_gates.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(aggregate12(), indent=1, default=str))
    print(json.dumps(aggregate11(), indent=1, default=str))
