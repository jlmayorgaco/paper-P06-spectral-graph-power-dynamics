# ruff: noqa: E501
"""E13 reduced / low-dimensional models and E14 decision certificate (prereg E13, E14)."""

from __future__ import annotations

import json
import time

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import E11_spectral as SP

PHASE13 = "E13"
PHASE14 = "E14"
HOLD = [f"H{i:02d}" for i in range(1, 25)]
GM_R = (2, 4, 8, 12, 16, 24, 39)
POD_R = (2, 4, 8, 12, 16)
TS = {"TSa": ("i_d", "i_q", "x_id", "x_iq"), "TSb": ("i_d", "i_q", "x_id", "x_iq", "p_filt", "q_filt"),
      "TSc": ("i_d", "i_q", "x_id", "x_iq", "p_filt", "q_filt", "theta_pll", "x_pll")}
CONTOUR = {"s0": 0.002, "s1": 20.0, "w1": 2 * np.pi * 5, "n": 400}
CERT_MAX = 0.9


def census_sample():
    rng = np.random.default_rng(20260923)
    subs = C.subsets(C.V9)
    idx = rng.choice(len(subs), size=64, replace=False)
    return [list(subs[i]) for i in sorted(idx)]


def pod_basis_path():
    return I.RESULTS / "CDW_E13_pod_basis.npy"


def build_pod_basis():
    """SVD of discovery critical-mode voltage shapes (D01-D15, core 16) from the E1 raw data."""

    cols = []
    for rec in I.Store("E01").all():
        t = rec["task"]
        if not rec.get("ok") or not t["pid"].startswith("D"):
            continue
        for r in rec["records"]:
            m = C.label(C.V4)
            del m
            mem = () if r["S"] == "BASE" else tuple(int(b) for b in r["S"].split("+"))
            if not set(mem) <= set(C.V4):
                continue
            for md in r.get("modes", []):
                if md["crit"]:
                    ph = np.asarray(md["phi"])
                    n = ph.size // 2
                    z = ph[:n] + 1j * ph[n:]
                    cols += [z.real, z.imag]
    X = np.array(cols).T  # 78 x K
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    np.save(pod_basis_path(), U[:, : max(POD_R)])
    return {"n_snapshots": X.shape[1], "sv_top": s[:10].tolist()}


def tasks13():
    out = [{"phase": PHASE13, "pid": p, "S": list(s), "set": "core"} for p in HOLD for s in C.subsets(C.V4)]
    out += [{"phase": PHASE13, "pid": p, "S": s, "set": "census"} for p in HOLD[:6] for s in census_sample()]
    return out


def fast_idx(labels, fam):
    names = TS[fam]
    return [k for k, lab in enumerate(labels) if "_gfl" in lab and lab.split("_gfl")[0] in names]


def transverse_alpha(A, r_x, w):
    tr = C.transverse_operator(A, r_x, w)
    ev = np.linalg.eigvals(tr.a_perp)
    return float(ev.real.max()), ev


def run_task13(task):
    theta = C.policy(task["pid"])
    t0 = time.perf_counter()
    try:
        case = C.solve(tuple(task["S"]), theta)
    except (C.InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as e:
        return {"status": "INFEASIBLE", "error": str(e)[:100]}
    a, d, jac = C.physical_matrices(case)
    r_x, _ = C.rotation_generator(case.dae, case.equilibrium.z)
    w = C.frequency_partner(case.dae).w
    te = time.perf_counter()
    alpha_full, _ = transverse_alpha(a, r_x, w)
    t_eig_full = time.perf_counter() - te
    status = C.classify_spectrum(C.transverse_operator(a, r_x, w).a_perp, C.transverse_operator(a, r_x, w).z.T @ d @ C.transverse_operator(a, r_x, w).z, C.SAFETY).status
    t_full = time.perf_counter() - t0
    out = {"status": status, "alpha_full": alpha_full, "n_x": int(a.shape[0]), "t_full": t_full, "t_eig_full": t_eig_full, "fam": {}}
    wg, U, _ = SP.graph_basis()
    Upod = np.load(pod_basis_path()) if pod_basis_path().exists() else None
    for fam, rs in (("GM", GM_R), ("POD", POD_R)):
        for r in rs:
            if fam == "GM":
                P = np.kron(U[:, :r], np.eye(2))
            else:
                if Upod is None:
                    continue
                P = Upod[:, :r]
            te = time.perf_counter()
            gzr = P.T @ jac.gz @ P
            Ar = jac.fx - jac.fz @ P @ np.linalg.solve(gzr, P.T @ jac.gx)
            al, _ = transverse_alpha(Ar, r_x, w)
            out["fam"][f"{fam}{r}"] = {"alpha": al, "state_ratio": 1.0, "net_dim_ratio": P.shape[1] / 78, "t_eig": time.perf_counter() - te}
    labels = case.dae.labels
    for fam in TS:
        f = fast_idx(labels, fam)
        if not f:
            out["fam"][fam] = {"alpha": alpha_full, "state_ratio": 1.0, "t_eig": t_eig_full, "trivial": True}
            continue
        s = [k for k in range(a.shape[0]) if k not in f]
        te = time.perf_counter()
        try:
            Ared = a[np.ix_(s, s)] - a[np.ix_(s, f)] @ np.linalg.solve(a[np.ix_(f, f)], a[np.ix_(f, s)])
            al, _ = transverse_alpha(Ared, r_x[s], w[s])
            out["fam"][fam] = {"alpha": al, "state_ratio": len(s) / a.shape[0], "t_eig": time.perf_counter() - te}
        except np.linalg.LinAlgError as e:
            out["fam"][fam] = {"alpha": None, "error": str(e)[:80]}
    return out


def aggregate13():
    rows = []
    for r in I.Store(PHASE13).all():
        if not r.get("ok") or r.get("status") == "INFEASIBLE":
            continue
        t = r["task"]
        for fam, v in r["fam"].items():
            rows.append({"pid": t["pid"], "set": t["set"], "S": C.label(t["S"]), "status_full": r["status"], "alpha_full": r["alpha_full"],
                         "family": fam, "alpha_red": v.get("alpha"), "state_ratio": v.get("state_ratio"), "net_dim_ratio": v.get("net_dim_ratio"),
                         "t_full": r["t_full"], "t_eig_full": r["t_eig_full"], "t_eig_red": v.get("t_eig"), "n_x": r["n_x"]})
    df = pd.DataFrame(rows)
    df["verdict_full"] = np.where(df.status_full == "UNSTABLE", 1, np.where(df.status_full == "STABLE", 0, -1))
    df["verdict_red"] = (df.alpha_red > 0).astype(int)
    df.to_csv(I.RESULTS / "CDW_E13_reduction_cases.csv", index=False)
    summ = []
    for fam, g in df.groupby("family"):
        g2 = g[g.verdict_full >= 0]
        acc = float((g2.verdict_full == g2.verdict_red).mean())
        # H exact per policy (core)
        hmatch = []
        taus = []
        for pid, gg in g[(g.set == "core")].groupby("pid"):
            full = {r.S: r.verdict_full for r in gg.itertuples()}
            red = {r.S: r.verdict_red for r in gg.itertuples()}
            hmatch.append(_H(full) == _H(red))
            a_f = {r.S: r.alpha_full for r in gg.itertuples()}
            a_r = {r.S: r.alpha_red for r in gg.itertuples()}
            df_, dr = [], []
            for S in C.subsets(C.V4):
                for i in C.V4:
                    if i in S:
                        continue
                    k1, k2 = C.label(S), C.label(S + (i,))
                    if k1 in a_f and k2 in a_f and a_r.get(k1) is not None and a_r.get(k2) is not None:
                        df_.append(a_f[k2] - a_f[k1])
                        dr.append(a_r[k2] - a_r[k1])
            taus.append(AN.kendall(df_, dr))
        summ.append({"family": fam, "verdict_acc": acc, "n": int(len(g2)), "H_exact_frac": float(np.mean(hmatch)) if hmatch else np.nan,
                     "kendall_marginals_median": float(np.nanmedian(taus)) if taus else np.nan,
                     "alpha_mae": float((g.alpha_red - g.alpha_full).abs().median()), "state_ratio_median": float(g.state_ratio.median()),
                     "eig_speedup_median": float((g.t_eig_full / g.t_eig_red).median()),
                     "end_to_end_speedup_median": float((g.t_full / (g.t_full - g.t_eig_full + g.t_eig_red)).median())})
    s = pd.DataFrame(summ)
    s.to_csv(I.RESULTS / "CDW_E13_reduction.csv", index=False)
    return s


def _H(verdicts):
    unsafe = [S for S, v in verdicts.items() if v == 1]
    mem = {S: set(() if S == "BASE" else map(int, S.split("+"))) for S in unsafe}
    minimal = [S for S in unsafe if not any(mem[r] < mem[S] for r in unsafe)]
    return "|".join(sorted(minimal)) or "EMPTY"


# ------------------------------------------------------------------ E14 --
def tasks14():
    return [{"phase": PHASE14, "pid": p, "S": list(s), "fam": fam} for p in HOLD for s in C.subsets(C.V4) if s for fam in TS]


def reduced_blocks(jac, labels, fam):
    f = fast_idx(labels, fam)
    s = [k for k in range(jac.fx.shape[0]) if k not in f]
    Ff = jac.fx[np.ix_(f, f)]
    fx = jac.fx[np.ix_(s, s)] - jac.fx[np.ix_(s, f)] @ np.linalg.solve(Ff, jac.fx[np.ix_(f, s)])
    fz = jac.fz[s] - jac.fx[np.ix_(s, f)] @ np.linalg.solve(Ff, jac.fz[f])
    gx = jac.gx[:, s] - jac.gx[:, f] @ np.linalg.solve(Ff, jac.fx[np.ix_(f, s)])
    gz = jac.gz - jac.gx[:, f] @ np.linalg.solve(Ff, jac.fz[f])
    return fx, fz, gx, gz, s


def contour():
    s0, s1, w1, n = CONTOUR["s0"], CONTOUR["s1"], CONTOUR["w1"], CONTOUR["n"]
    return np.concatenate([np.linspace(s0, s1, n) - 1j * w1, s1 + 1j * np.linspace(-w1, w1, n),
                           np.linspace(s1, s0, n) + 1j * w1, s0 + 1j * np.linspace(w1, -w1, n)])


def inside(ev):
    return (ev.real > CONTOUR["s0"]) & (ev.real < CONTOUR["s1"]) & (np.abs(ev.imag) < CONTOUR["w1"])


def run_task14(task):
    theta = C.policy(task["pid"])
    case = C.solve(tuple(task["S"]), theta)
    a, d, jac = C.physical_matrices(case)
    labels = case.dae.labels
    r_x, _ = C.rotation_generator(case.dae, case.equilibrium.z)
    w = C.frequency_partner(case.dae).w
    fx, fz, gx, gz, s = reduced_blocks(jac, labels, task["fam"])
    ev_full = np.linalg.eigvals(C.transverse_operator(a, r_x, w).a_perp)
    Ared = fx - fz @ np.linalg.solve(gz, gx)
    ev_red = np.linalg.eigvals(C.transverse_operator(Ared, r_x[s], w[s]).a_perp)
    n_full, n_red = int(inside(ev_full).sum()), int(inside(ev_red).sum())
    rhp_outside = int(((ev_full.real > 0) & ~inside(ev_full) & (np.abs(ev_full) > 1e-6)).sum())
    poles = np.concatenate([np.linalg.eigvals(jac.fx), np.linalg.eigvals(fx)])
    n_pole_in = int(inside(poles).sum())
    out = {"n_full": n_full, "n_red": n_red, "verdict_full": int(n_full > 0), "verdict_red": int(n_red > 0),
           "n_pole_inside": n_pole_in, "rhp_outside": rhp_outside}
    if n_pole_in or rhp_outside:
        out["cert"] = "ABSTAIN"
        out["reason"] = "pole_inside" if n_pole_in else "scope"
        return out
    nx_, nr = jac.fx.shape[0], fx.shape[0]
    worst = 0.0
    t0 = time.perf_counter()
    for z in contour():
        T = jac.gz + jac.gx @ np.linalg.solve(z * np.eye(nx_) - jac.fx, jac.fz)
        Tr = gz + gx @ np.linalg.solve(z * np.eye(nr) - fx, fz)
        M = np.linalg.solve(Tr, T - Tr)
        worst = max(worst, float(np.linalg.norm(M, 2)))
    out.update(cert_value=worst, cert="SAMPLED-CERTIFIED" if worst < CERT_MAX else "ABSTAIN",
               reason="" if worst < CERT_MAX else "bound", t_contour=time.perf_counter() - t0)
    return out


def aggregate14():
    rows = []
    for r in I.Store(PHASE14).all():
        if not r.get("ok"):
            continue
        t = r["task"]
        rows.append({"pid": t["pid"], "S": C.label(t["S"]), "fam": t["fam"], **{k: v for k, v in r.items() if k not in ("task", "ok")}})
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E14_certificate.csv", index=False)
    summ = []
    for fam, g in df.groupby("fam"):
        cert = g[g.cert == "SAMPLED-CERTIFIED"]
        summ.append({"fam": fam, "n": int(len(g)), "coverage": float(len(cert) / len(g)), "false_certifications": int((cert.verdict_full != cert.verdict_red).sum()),
                     "abstain_pole": int((g.reason == "pole_inside").sum()), "abstain_scope": int((g.reason == "scope").sum()),
                     "abstain_bound": int((g.reason == "bound").sum()), "reduced_verdict_acc_all": float((g.verdict_full == g.verdict_red).mean()),
                     "reduced_verdict_acc_certified": float((cert.verdict_full == cert.verdict_red).mean()) if len(cert) else np.nan})
    s = pd.DataFrame(summ)
    s.to_csv(I.RESULTS / "CDW_E14_summary.csv", index=False)
    return s


def gold_c():
    r13 = pd.read_csv(I.RESULTS / "CDW_E13_reduction.csv")
    r14 = pd.read_csv(I.RESULTS / "CDW_E14_summary.csv") if (I.RESULTS / "CDW_E14_summary.csv").exists() else pd.DataFrame()
    res = []
    for row in r13.itertuples():
        cert = r14[r14.fam == row.family] if len(r14) else pd.DataFrame()
        abst = float(1 - cert.coverage.iloc[0]) if len(cert) else 0.0
        acc = float(cert.reduced_verdict_acc_certified.iloc[0]) if len(cert) and np.isfinite(cert.reduced_verdict_acc_certified.iloc[0]) else row.verdict_acc
        false_c = int(cert.false_certifications.iloc[0]) if len(cert) else 0
        ok = (row.state_ratio_median <= 0.70 and row.end_to_end_speedup_median >= 3 and acc >= 0.99 and abst <= 0.20
              and false_c == 0 and row.H_exact_frac >= 0.90)
        res.append({"family": row.family, "state_ratio": row.state_ratio_median, "e2e_speedup": row.end_to_end_speedup_median,
                    "eig_speedup": row.eig_speedup_median, "acc": acc, "abstention": abst, "false_cert": false_c,
                    "H_exact": row.H_exact_frac, "pass": bool(ok)})
    g = {"families": res, "GOLD_C_pass": bool(any(r["pass"] for r in res))}
    I.atomic_write_json(I.RESULTS / "CDW_E13_E14_gates.json", g)
    return g


if __name__ == "__main__":
    print(aggregate13())
    print(aggregate14())
    print(json.dumps(gold_c(), indent=1, default=str))
