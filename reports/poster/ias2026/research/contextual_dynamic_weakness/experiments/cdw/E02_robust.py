# ruff: noqa: E501
"""E2 robust contextuality under the TX4 (discovery) and CDW (holdout) envelope draws."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

PHASE = "E02"
ENVS = ("EM-f", "EM-u", "EC", "EMC")
CORE16 = C.subsets(C.V4)


def verify_tx4_draws() -> dict:
    """Prereg rule: regenerated TX4 draws must reproduce the frozen PCV05 factors exactly."""

    frozen = pd.read_csv(C.RESEARCH / "results/PCV/PCV05/PCV05_draws.csv")
    frozen = frozen[frozen.envelope.isin(ENVS)]
    reg = {(d["envelope"], d["draw"]): d for d in C.tx4_draws()}
    n = bad = 0
    for r in frozen.itertuples():
        f = json.loads(r.factors)
        d = reg.get((r.envelope, int(r.draw)))
        n += 1
        if d is None:
            bad += 1
            continue
        same = (f["fleet"] == d["fleet"] and {str(k): v for k, v in f["unit"].items()} == d["unit"]
                and {str(k): v for k, v in f["conv"].items()} == d["conv"])
        bad += not same
    return {"n_checked": n, "n_mismatch": bad, "ok": bad == 0 and n == 400}


def tasks():
    out = []
    ver = verify_tx4_draws()
    I.atomic_write_json(I.RESULTS / "CDW_E2_tx4_draw_verification.json", ver)
    tx4 = C.tx4_draws() if ver["ok"] else []
    cdw = C.cdw_draws()
    for pid in ("D01", "D03"):
        out.append({"phase": PHASE, "kind": "core", "pid": pid, "source": "NOMINAL", "env": "NOMINAL", "draw": 0})
        for d in tx4:
            out.append({"phase": PHASE, "kind": "core", "pid": pid, "source": "TX4", "env": d["envelope"], "draw": d["draw"]})
        for d in cdw:
            out.append({"phase": PHASE, "kind": "core", "pid": pid, "source": "CDW", "env": d["envelope"], "draw": d["draw"]})
    subs = C.subsets(C.V9)
    for env in ENVS:
        for k in range(10):
            for c in range(0, 512, 64):
                out.append({"phase": PHASE, "kind": "census", "pid": "D01", "source": "CDW", "env": env, "draw": k,
                            "subsets": [list(s) for s in subs[c: c + 64]]})
    return out


def _draw(task):
    if task["source"] == "NOMINAL":
        return None
    pool = C.tx4_draws() if task["source"] == "TX4" else C.cdw_draws()
    for d in pool:
        if d["envelope"] == task["env"] and d["draw"] == task["draw"]:
            return d
    raise KeyError(task)


def run_task(task):
    theta = C.policy(task["pid"])
    draw = _draw(task)
    subs = CORE16 if task["kind"] == "core" else [tuple(s) for s in task["subsets"]]
    recs = []
    for s in subs:
        r = C.eval_portfolio(s, theta, modes=False, draw=draw)
        recs.append({"S": C.label(s), **{k: r.get(k) for k in ("status", "alpha", "lam_hz", "rhp")}})
    return {"records": recs}


def lattice_metrics(alpha: dict, status: dict, cands=C.V4, tau=C.TAU_MAT):
    """Hypergraph, kappa and marginal sign structure of one lattice."""

    st = {s: status[C.label(s)] for s in C.subsets(cands)}
    base_ok = st[()] == "STABLE"
    unsafe = [s for s, v in st.items() if v == "UNSTABLE"]
    minimal = sorted((s for s in unsafe if not any(set(r) < set(s) for r in unsafe)), key=lambda s: (len(s), s))
    H = "|".join(C.label(e) for e in minimal) or "EMPTY"
    marg = {}
    for S in C.subsets(cands):
        for i in cands:
            if i in S:
                continue
            a, b = alpha[C.label(S)], alpha[C.label(S + (i,))]
            marg[(C.label(S), i)] = (b - a, st[S] == "STABLE")
    rev_stable = False
    for i in cands:
        cls = [AN.sign_class(d, tau) for (S, j), (d, stab) in marg.items() if j == i and stab]
        if -1 in cls and 1 in cls:
            rev_stable = True
    return {"H": H if base_ok else "BASE_UNSTABLE", "kappa": min((len(e) for e in minimal), default=np.inf) if base_ok else np.nan,
            "rev_stable": rev_stable, "marg": marg}


def aggregate():
    recs = [r for r in I.Store(PHASE).all()]
    bad = [r for r in recs if not r.get("ok")]
    if bad:
        raise RuntimeError(f"{len(bad)} failed E02 tasks: {bad[0].get('error')}")
    core = [r for r in recs if r["task"]["kind"] == "core"]
    rows = []
    nominal = {}
    for r in core:
        t = r["task"]
        alpha = {x["S"]: x["alpha"] for x in r["records"]}
        status = {x["S"]: x["status"] for x in r["records"]}
        if any(v in ("INFEASIBLE",) for v in status.values()):
            rows.append({**{k: t[k] for k in ("pid", "source", "env", "draw")}, "feasible": False})
            continue
        m = lattice_metrics(alpha, status)
        rec = {**{k: t[k] for k in ("pid", "source", "env", "draw")}, "feasible": True, "H": m["H"], "kappa": m["kappa"],
               "rev_stable": m["rev_stable"], "alpha_H4": alpha["30+33+35+37"]}
        rec["_marg"] = m["marg"]
        rows.append(rec)
        if t["source"] == "NOMINAL":
            nominal[t["pid"]] = m
    df = pd.DataFrame(rows)
    out = []
    for r in df.itertuples():
        d = r._asdict()
        if not d.get("feasible") or r.pid not in nominal:
            out.append({k: v for k, v in d.items() if k != "_marg"})
            continue
        nm = nominal[r.pid]
        marg = d["_marg"]
        keys = [k for k in nm["marg"] if k in marg]
        same_sign = [AN.sign_class(marg[k][0], C.TAU_MAT) == AN.sign_class(nm["marg"][k][0], C.TAU_MAT) for k in keys]
        # per-context ranking persistence (Kendall over available interventions, contexts with >= 3)
        taus, top1 = [], []
        for S in {k[0] for k in keys}:
            av = [k for k in keys if k[0] == S]
            if len(av) >= 3:
                a = [marg[k][0] for k in av]
                b = [nm["marg"][k][0] for k in av]
                taus.append(AN.kendall(a, b))
                top1.append(int(np.argmin(a) == np.argmin(b)))
        d.update(sign_persistence=float(np.mean(same_sign)), kendall_ctx=float(np.nanmean(taus)) if taus else np.nan,
                 top1_ctx=float(np.mean(top1)) if top1 else np.nan, witness_changed=r.H != nm["H"])
        out.append({k: v for k, v in d.items() if k != "_marg"})
    res = pd.DataFrame(out).drop(columns=["Index"], errors="ignore")
    res.to_parquet(I.RESULTS / "CDW_E2_robust_contextuality.parquet", index=False)
    summ = (res[res.source != "NOMINAL"].groupby(["pid", "source", "env"])
            .agg(n=("feasible", "size"), n_feasible=("feasible", "sum"), cov_rev_stable=("rev_stable", "mean"),
                 cov_witness_changed=("witness_changed", "mean"), sign_persistence=("sign_persistence", "mean"),
                 kendall_ctx=("kendall_ctx", "mean"), top1_ctx=("top1_ctx", "mean"))
            .reset_index())
    cond = (res[(res.source != "NOMINAL") & (res.witness_changed == True)]  # noqa: E712
            .groupby(["pid", "source", "env"]).rev_stable.mean().rename("cov_rev_given_witness_changed").reset_index())
    summ = summ.merge(cond, on=["pid", "source", "env"], how="left")
    summ.to_csv(I.RESULTS / "CDW_E2_robust_summary.csv", index=False)
    a4 = summ[(summ.pid == "D01") & (summ.source == "CDW")]
    gate = {"A4_cov_by_env": dict(zip(a4.env, a4.cov_rev_stable.round(4), strict=True)),
            "A4_envs_ge_0.5": int((a4.cov_rev_stable >= 0.5).sum())}
    gate["A4_pass"] = bool(gate["A4_envs_ge_0.5"] >= 3)
    # census on 10 CDW draws per envelope (descriptive): reversal presence on V9
    cen = [r for r in recs if r["task"]["kind"] == "census"]
    by = {}
    for r in cen:
        t = r["task"]
        by.setdefault((t["env"], t["draw"]), []).extend(r["records"])
    crow = []
    for (env, k), lst in by.items():
        alpha = {x["S"]: x["alpha"] for x in lst}
        status = {x["S"]: x["status"] for x in lst}
        if len(alpha) < 512 or any(s == "INFEASIBLE" for s in status.values()):
            crow.append({"env": env, "draw": k, "complete": False})
            continue
        n_rev = 0
        for i in C.V9:
            cls = []
            for S in C.subsets(C.V9):
                if i in S or status[C.label(S)] != "STABLE":
                    continue
                cls.append(AN.sign_class(alpha[C.label(S + (i,))] - alpha[C.label(S)], C.TAU_MAT))
            n_rev += (-1 in cls) and (1 in cls)
        crow.append({"env": env, "draw": k, "complete": True, "n_reversing_interventions": n_rev, "rev_stable": n_rev > 0})
    cdf = pd.DataFrame(crow)
    cdf.to_csv(I.RESULTS / "CDW_E2_census_draws.csv", index=False)
    gate["census_cov_by_env"] = cdf[cdf.complete].groupby("env").rev_stable.mean().round(4).to_dict() if len(cdf) else {}
    I.atomic_write_json(I.RESULTS / "CDW_E2_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
