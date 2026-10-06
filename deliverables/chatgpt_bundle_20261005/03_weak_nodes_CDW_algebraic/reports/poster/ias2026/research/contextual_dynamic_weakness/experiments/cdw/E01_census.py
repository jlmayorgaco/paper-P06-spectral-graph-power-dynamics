# ruff: noqa: E501
"""E1 contextual sign-reversal census (+ E1b submodularity, E10 Shapley aggregates).

Prereg: docs/CDW_PREREG_V1.md section 2 (E1, E1b, E10).
"""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

PHASE = "E01"
CHUNK = 32
POLICIES = list(C.DISCOVERY) + [f"H{i:02d}" for i in range(1, 25)]


def tasks(policies=POLICIES, device="gfl"):
    subs = C.subsets(C.V9)
    out = []
    for pid in policies:
        for c in range(0, len(subs), CHUNK):
            out.append({"phase": PHASE, "pid": pid, "device": device, "subsets": [list(s) for s in subs[c: c + CHUNK]]})
    return out


def run_task(task):
    theta = C.policy(task["pid"])
    recs = []
    for s in task["subsets"]:
        r = C.eval_portfolio(tuple(s), theta, modes=True, device=task["device"])
        r["S"] = C.label(s)
        recs.append(r)
    return {"records": recs}


# ------------------------------------------------------------------ aggregate --
def load(phase=PHASE, device="gfl"):
    rows, modes = [], {}
    for rec in I.Store(phase).all():
        if not rec.get("ok"):
            raise RuntimeError(f"failed task in {phase}: {rec.get('error')}")
        t = rec["task"]
        if t.get("device", "gfl") != device:
            continue
        for r in rec["records"]:
            key = (t["pid"], r["S"])
            modes[key] = [compact(m) for m in r.get("modes", [])]
            rows.append({"pid": t["pid"], "S": r["S"], "size": 0 if r["S"] == "BASE" else r["S"].count("+") + 1,
                         **{k: r.get(k) for k in ("status", "alpha", "lam_re", "lam_hz", "rhp", "gap2", "n_x", "error")}})
    df = pd.DataFrame(rows)
    return df, modes


def compact(m):
    ph = np.asarray(m["phi"], np.float32)
    n = ph.size // 2
    return {"re": m["re"], "hz": m["hz"], "crit": m["crit"], "phi": (ph[:n] + 1j * ph[n:]).astype(np.complex64)}


def mac_c(a, b) -> float:
    num = abs(np.vdot(a, b)) ** 2
    den = float(np.vdot(a, a).real * np.vdot(b, b).real) or 1e-300
    return float(num / den)


def members_of(label):
    return () if label == "BASE" else tuple(int(b) for b in label.split("+"))


def crit_mode(ms):
    for m in ms:
        if m["crit"]:
            return m
    return None


def match(mode, cands):
    """Best MAC match of `mode` among cands within DF_MAX; returns (mode, mac) or (None, 0)."""

    best, bm = None, 0.0
    for c in cands:
        if abs(c["hz"] - mode["hz"]) > C.DF_MAX:
            continue
        m = mac_c(mode["phi"], c["phi"])
        if m > bm:
            best, bm = c, m
    if bm < C.MAC_MIN:
        return None, bm
    return best, bm


def marginals(df, modes, cands=C.V9):
    out = []
    alpha = {(r.pid, r.S): r.alpha for r in df.itertuples()}
    status = {(r.pid, r.S): r.status for r in df.itertuples()}
    for pid in sorted(df.pid.unique()):
        for S in C.subsets(cands):
            lS = C.label(S)
            aS = alpha.get((pid, lS))
            for i in cands:
                if i in S:
                    continue
                lSi = C.label(S + (i,))
                aSi = alpha.get((pid, lSi))
                rec = {"pid": pid, "S": lS, "size": len(S), "i": i, "alpha_S": aS, "alpha_Si": aSi,
                       "delta": (aSi - aS) if (aS is not None and aSi is not None) else np.nan,
                       "stable_S": status.get((pid, lS)) == "STABLE",
                       "status_S": status.get((pid, lS)), "status_Si": status.get((pid, lSi))}
                cm = crit_mode(modes.get((pid, lS), []))
                cmi = crit_mode(modes.get((pid, lSi), []))
                if cm is not None and modes.get((pid, lSi)):
                    m, mc = match(cm, modes[(pid, lSi)])
                    rec["mac"] = mc
                    if m is not None:
                        rec["tracked_delta"] = m["re"] - cm["re"]
                        rec["same_mode"] = bool(m["crit"])
                    else:
                        rec["tracked_delta"] = np.nan
                        rec["same_mode"] = False
                else:
                    rec["mac"], rec["tracked_delta"], rec["same_mode"] = np.nan, np.nan, False
                rec["crit_hz_S"] = cm["hz"] if cm else np.nan
                rec["crit_hz_Si"] = cmi["hz"] if cmi else np.nan
                out.append(rec)
    return pd.DataFrame(out)


def per_intervention(mg, tau=C.TAU_MAT):
    rows = []
    for (pid, i), g in mg.groupby(["pid", "i"]):
        d = g.delta.to_numpy(float)
        ok = np.isfinite(d)
        d = d[ok]
        gs = g[ok]
        cls = np.array([AN.sign_class(x, tau) for x in d])
        st = gs.stable_S.to_numpy(bool)
        same = gs.same_mode.to_numpy(bool)
        counts = [(cls == -1).sum(), (cls == 0).sum(), (cls == 1).sum()]

        def rev(mask):
            return bool(((cls == -1) & mask).any() and ((cls == 1) & mask).any())

        rec = {"pid": pid, "i": i, "n": int(d.size), "frac_stab": counts[0] / max(d.size, 1), "frac_neutral": counts[1] / max(d.size, 1),
               "frac_destab": counts[2] / max(d.size, 1), "min": float(d.min()) if d.size else np.nan,
               "median": float(np.median(d)) if d.size else np.nan, "max": float(d.max()) if d.size else np.nan,
               "var": float(d.var()) if d.size else np.nan, "entropy": AN.entropy3(counts),
               "n_stab": int(counts[0]), "n_destab": int(counts[2]),
               "rev_global": rev(np.ones_like(st)), "rev_stable": rev(st), "rev_tracked": rev(st & same)}
        # strongest stable-context reversal pair
        if rec["rev_stable"]:
            neg = gs[(cls == -1) & st]
            pos = gs[(cls == 1) & st]
            a, b = neg.loc[neg.delta.idxmin()], pos.loc[pos.delta.idxmax()]
            rec.update(strongest=min(abs(a.delta), abs(b.delta)), S_stab=a.S, S_destab=b.S,
                       d_stab=float(a.delta), d_destab=float(b.delta))
        # strict variant: same family across the two contexts
        rec["rev_tracked_samefamily"] = False
        rows.append(rec)
    return pd.DataFrame(rows)


def node_only_test(mg, df, cands=C.V9, tau=C.TAU_MAT):
    """Optimal fixed ranking (DP) on the policy's own stable contexts; agreement p* and regret."""

    rows = []
    idx = {b: k for k, b in enumerate(cands)}
    for pid, g in mg.groupby("pid"):
        g = g[g.stable_S & np.isfinite(g.delta)]
        W = np.zeros((len(cands), len(cands)))
        ctx = {}
        for S, gg in g.groupby("S"):
            ctx[S] = dict(zip(gg.i, gg.delta, strict=True))
        total = 0
        for d in ctx.values():
            av = list(d)
            for a in av:
                for b in av:
                    if a != b and d[a] < d[b] - tau:  # a better than b (more stabilizing)
                        W[idx[a], idx[b]] += 1
                        total += 1
        order, val = AN.optimal_linear_order(W)
        pstar = val / total if total else np.nan
        rank = {cands[j]: r for r, j in enumerate(order)}
        n_ctx = n_bad = 0
        regrets = []
        for d in ctx.values():
            if len(d) < 2:
                continue
            n_ctx += 1
            choice = min(d, key=lambda u: rank[u])
            best = min(d.values())
            reg = d[choice] - best
            regrets.append(reg)
            n_bad += reg >= tau
        rows.append({"pid": pid, "p_star": pstar, "n_triples": total, "order": "-".join(str(cands[j]) for j in order),
                     "n_ctx": n_ctx, "frac_regret": n_bad / n_ctx if n_ctx else np.nan,
                     "median_regret": float(np.median(regrets)) if regrets else np.nan,
                     "max_regret": float(np.max(regrets)) if regrets else np.nan})
    return pd.DataFrame(rows)


def aggregate():
    df, modes = load()
    base = df[df.S == "BASE"].set_index("pid").alpha
    df.to_parquet(I.RESULTS / "CDW_E1_portfolios.parquet", index=False)
    mg = marginals(df, modes)
    mg.to_parquet(I.RESULTS / "CDW_E1_contextual_marginals.parquet", index=False)
    pi = per_intervention(mg)
    pi["base_alpha"] = pi.pid.map(base)
    pi.to_csv(I.RESULTS / "CDW_E1_summary.csv", index=False)
    no = node_only_test(mg, df)
    no["base_alpha"] = no.pid.map(base)
    pol = pi.groupby("pid").agg(rev_global=("rev_global", "any"), rev_stable=("rev_stable", "any"),
                                rev_tracked=("rev_tracked", "any"), n_rev_interventions=("rev_stable", "sum")).reset_index()
    pol = pol.merge(no, on="pid")
    pol["base_stable"] = pol.base_alpha < 0
    pol["n_infeasible"] = pol.pid.map(df.groupby("pid").apply(lambda g: int((g.status == "INFEASIBLE").sum()), include_groups=False))
    pol["n_unresolved"] = pol.pid.map(df.groupby("pid").apply(lambda g: int((g.status == "BOUNDARY_OR_UNRESOLVED").sum()), include_groups=False))
    pol["split"] = np.where(pol.pid.str.startswith("H"), "holdout", "discovery")
    pol.to_csv(I.RESULTS / "CDW_E1_policy_summary.csv", index=False)
    # tau_res sensitivity
    pi_res = per_intervention(mg, tau=C.TAU_RES)
    pi_res.to_csv(I.RESULTS / "CDW_E1_summary_tau_res.csv", index=False)
    h = pol[(pol.split == "holdout") & pol.base_stable]
    gates = {
        "n_holdout_base_stable": int(len(h)), "n_holdout_base_unstable": int(((pol.split == "holdout") & ~pol.base_stable).sum()),
        "A1_frac": float(h.rev_stable.mean()) if len(h) else np.nan,
        "A2_frac": float(h.rev_tracked.mean()) if len(h) else np.nan,
        "A3_median_pstar": float(h.p_star.median()) if len(h) else np.nan,
        "A3_median_frac_regret": float(h.frac_regret.median()) if len(h) else np.nan,
    }
    gates["A1_pass"] = bool(gates["A1_frac"] >= 0.75)
    gates["A2_pass"] = bool(gates["A2_frac"] >= 0.50)
    gates["A3_pass"] = bool(gates["A3_median_pstar"] <= 0.90 and gates["A3_median_frac_regret"] >= 0.10)
    gates["H1_refuted_stop"] = bool(gates["A1_frac"] < 0.25)
    I.atomic_write_json(I.RESULTS / "CDW_E1_gates.json", gates)
    aggregate_E1b(df, modes)
    aggregate_E10(df)
    return gates


def aggregate_E1b(df, modes):
    rows, cex = [], []
    for pid, g in df.groupby("pid"):
        v = {members_of(r.S): r.alpha for r in g.itertuples() if np.isfinite(r.alpha)}
        sd = AN.second_differences(v, C.V9)
        d = np.array([x[3] for x in sd])
        sub_v = d > C.TAU_RES
        sup_v = d < -C.TAU_RES
        # tracked margin: mode matching the base critical mode
        base_c = crit_mode(modes.get((pid, "BASE"), []))
        vt = {}
        if base_c is not None:
            for r in g.itertuples():
                m, _ = match(base_c, modes.get((pid, r.S), []))
                if m is not None:
                    vt[members_of(r.S)] = m["re"]
        sdt = AN.second_differences(vt, C.V9)
        dt = np.array([x[3] for x in sdt]) if sdt else np.array([])
        rows.append({"pid": pid, "n_squares": int(d.size), "n_submod_viol": int(sub_v.sum()), "n_supermod_viol": int(sup_v.sum()),
                     "max_submod_viol": float(d.max()) if d.size else np.nan, "max_supermod_viol": float(-d.min()) if d.size else np.nan,
                     "submodular": bool(not sub_v.any()), "supermodular": bool(not sup_v.any()),
                     "tracked_n_squares": int(dt.size), "tracked_n_submod_viol": int((dt > C.TAU_RES).sum()),
                     "tracked_n_supermod_viol": int((dt < -C.TAU_RES).sum()),
                     "tracked_frac_viol": float(((dt > C.TAU_RES) | (dt < -C.TAU_RES)).mean()) if dt.size else np.nan,
                     "primary_frac_viol": float((sub_v | sup_v).mean()) if d.size else np.nan})
        for r in range(0, 8):
            sel = [x for x in sd if len(x[0]) == r]
            if not sel:
                continue
            a = max(sel, key=lambda x: x[3])
            b = min(sel, key=lambda x: x[3])
            for kind, x in (("submod_violation", a), ("supermod_violation", b)):
                if (kind == "submod_violation" and x[3] > C.TAU_RES) or (kind == "supermod_violation" and x[3] < -C.TAU_RES):
                    cex.append({"pid": pid, "kind": kind, "S": C.label(x[0]), "i": x[1], "j": x[2], "size_S": r, "second_diff": x[3]})
    pd.DataFrame(rows).to_csv(I.RESULTS / "CDW_E1b_submodularity.csv", index=False)
    pd.DataFrame(cex).to_csv(I.RESULTS / "CDW_E1b_minimal_counterexamples.csv", index=False)


def aggregate_E10(df):
    mg = pd.read_parquet(I.RESULTS / "CDW_E1_contextual_marginals.parquet")
    rows = []
    for pid, g in df.groupby("pid"):
        v = {members_of(r.S): r.alpha for r in g.itertuples()}
        if any(not np.isfinite(x) for x in v.values()):
            continue
        a0 = v[()]
        vv = {k: x - a0 for k, x in v.items()}
        phi = AN.shapley(vv, C.V9)
        for i in C.V9:
            d = mg[(mg.pid == pid) & (mg.i == i)].delta.to_numpy(float)
            d = d[np.isfinite(d)]
            cls = np.array([AN.sign_class(x, C.TAU_MAT) for x in d])
            rev = bool((cls == -1).any() and (cls == 1).any())
            rows.append({"pid": pid, "i": i, "shapley": phi[i], "mean_marginal": float(d.mean()), "std_marginal": float(d.std()),
                         "entropy": AN.entropy3([(cls == -1).sum(), (cls == 0).sum(), (cls == 1).sum()]), "reversal": rev,
                         "hidden_reversal": bool(abs(phi[i]) <= C.TAU_MAT and d.std() >= 3 * C.TAU_MAT and rev)})
    pd.DataFrame(rows).to_csv(I.RESULTS / "CDW_E10_shapley.csv", index=False)


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1))
