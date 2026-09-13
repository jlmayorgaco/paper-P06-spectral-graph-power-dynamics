# ruff: noqa: E501
"""H31 exploratory computations requested by the internal reviewers (POST-HOC; no prereg verdict changes).

(1) participation factors of the fast real critical eigenvalues (new holdout, base-stable policies);
(2) fractional-replacement sweeps rho_i in [0, 1] into a fast portfolio: alpha, sigma_min(g_z);
(3) continuation check of level-C/D tracked-mode identity: sweep rho_i in [0, 1] for each endpoint
    marginal and follow the critical mode of S continuously; compare with the MAC-matched mode.
(4) timing of the total derivative vs finite re-solve per branch at one condition."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json
import time

import numpy as np
import pandas as pd

import _cdw as C
import _infra as I
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

PHASE = "H_H31x"
R = HI.RESULTS
SEED = 20260934
SG_STATES = ["delta", "omega", "eq1", "ed1", "efd", "pw", "pl", "avr_ll"]
GFL_STATES = ["theta", "x_pll", "p_f", "q_f", "x_p", "x_q", "i_d", "i_q", "x_id", "x_iq", "x_v"]


def solve_frac(members, i, rho, theta):
    plan = ReplacementPlan.of({**{b: 1.0 for b in members}, i: rho})
    return solve_case(plan, **C.case_kwargs(theta))


def state_labels(case):
    lab = [None] * case.dae.n_x
    for s in case.dae.slots:
        names = SG_STATES if s.kind == "sg" else GFL_STATES
        for k in range(s.start, s.stop):
            lab[k] = (s.kind, s.bus, names[k - s.start] if k - s.start < len(names) else f"s{k - s.start}")
    return lab


def participation(case):
    a, _, jac = C.physical_matrices(case)
    vals, vl, vr = __import__("scipy").linalg.eig(a, left=True, right=True)
    ok = np.abs(vals) > C.STRUCT_ZERO
    j = np.where(ok)[0][np.argmax(vals[ok].real)]
    p = np.abs(vl[:, j].conj() * vr[:, j])
    p = p / p.sum()
    lab = state_labels(case)
    groups = {}
    for k, w in enumerate(p):
        kind, bus, nm = lab[k]
        groups[f"{kind}:{nm}"] = groups.get(f"{kind}:{nm}", 0.0) + float(w)
    sv = np.linalg.svd(jac.gz, compute_uv=False)
    return {"lam": [float(vals[j].real), float(vals[j].imag)], "groups": dict(sorted(groups.items(), key=lambda kv: -kv[1])[:8]),
            "sigma_min_gz": float(sv.min()), "cond_gz": float(sv.max() / sv.min())}


def tasks():
    mg = pd.read_parquet(R / "H03_marginals.parquet", columns=["pid", "mask", "i", "stable_S", "alpha_S", "alpha_Si", "status_Si", "hz_Si", "lvA", "lvC"])
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    stable_new = set(pol[(pol.split == "new") & pol.base_stable].pid)
    rng = np.random.default_rng(SEED)
    out = []
    fast = mg[mg.pid.isin(stable_new) & mg.lvA & (mg.alpha_Si >= 10)]
    pick = fast.iloc[rng.choice(len(fast), size=30, replace=False)]
    for r in pick.itertuples():
        out.append({"phase": PHASE, "kind": "part", "pid": r.pid, "mask": int(r.mask), "i": int(r.i)})
    for r in pick.iloc[:6].itertuples():
        out.append({"phase": PHASE, "kind": "sweep", "pid": r.pid, "mask": int(r.mask), "i": int(r.i)})
    pairs = pd.read_parquet(R / "H03_nested_pairs.parquet")
    pc = pairs[(pairs.level == "C") & pairs.pid.isin(stable_new)]
    rev = json.loads((R / "H31_revision.json").read_text())
    ex = pd.DataFrame(rev["levelD_examples"])
    samp = pc.iloc[rng.choice(len(pc), size=30, replace=False)]
    for r in pd.concat([ex[["pid", "i", "S1", "S2"]], samp[["pid", "i", "S1", "S2"]]]).itertuples():
        for S in (r.S1, r.S2):
            out.append({"phase": PHASE, "kind": "cont", "pid": r.pid, "S": S, "i": int(r.i)})
    out.append({"phase": PHASE, "kind": "timing", "pid": "HARDENING_H01"})
    return out


def members_of(label):
    return () if label == "BASE" else tuple(int(u) for u in label.split("+"))


def run_task(t):
    theta = HI.theta_of(t["pid"])
    if t["kind"] == "part":
        S = tuple(C.V9[b] for b in range(9) if t["mask"] >> b & 1)
        case = C.solve(S + (t["i"],), theta)
        return participation(case)
    if t["kind"] == "sweep":
        S = tuple(C.V9[b] for b in range(9) if t["mask"] >> b & 1)
        rows = []
        for rho in np.linspace(0.0, 1.0, 41):
            try:
                case = solve_frac(S, t["i"], float(rho), theta)
            except Exception as e:  # noqa: BLE001
                rows.append({"rho": float(rho), "error": str(e)[:80]})
                continue
            a, _, jac = C.physical_matrices(case)
            ev = np.linalg.eigvals(a)
            ev = ev[np.abs(ev) > C.STRUCT_ZERO]
            sv = np.linalg.svd(jac.gz, compute_uv=False)
            real = ev[np.abs(ev.imag) < 1e-9]
            rows.append({"rho": float(rho), "alpha": float(ev.real.max()), "max_real_eig": float(real.real.max()) if real.size else None,
                         "min_real_eig": float(real.real.min()) if real.size else None, "sigma_min_gz": float(sv.min())})
        return {"rows": rows}
    if t["kind"] == "cont":
        S = members_of(t["S"])
        rec0 = C.eval_portfolio(S, theta, modes=False)
        lam0 = complex(rec0["lam_re"], 2 * np.pi * rec0["lam_hz"])
        track = lam0
        path = []
        for rho in np.linspace(0.0, 1.0, 21)[1:]:
            case = solve_frac(S, t["i"], float(rho), theta)
            ev = np.linalg.eigvals(C.physical_matrices(case)[0])
            ev = ev[ev.imag >= -1e-12]
            k = int(np.argmin(np.abs(ev - track)))
            gap_next = float(np.sort(np.abs(ev - track))[1]) if ev.size > 1 else np.inf
            track = complex(ev[k])
            path.append([float(rho), track.real, track.imag, abs(ev[k] - track), gap_next])
        rec1 = C.eval_portfolio(S + (t["i"],), theta, modes=False)
        crit1 = complex(rec1["lam_re"], 2 * np.pi * rec1["lam_hz"])
        return {"lam0": [lam0.real, lam0.imag], "tracked_end": [track.real, track.imag], "crit1": [crit1.real, crit1.imag],
                "continuation_matches_crit": bool(abs(track - crit1) < 1e-6 or abs(track - crit1.conjugate()) < 1e-6), "path": path}
    if t["kind"] == "timing":
        import _sens as SN
        import E34_sens as E34

        eng = SN.Engine(C.V4, theta)
        specs = E34.link_specs()[:6]
        t0 = time.perf_counter()
        eng.derivatives(specs, "SPR")
        t_der = (time.perf_counter() - t0) / len(specs)
        t0 = time.perf_counter()
        for sp in specs:
            SN.spr_case(eng, sp, 1.5)
            C.transverse_eval(SN.spr_case(eng, sp, 1.5), modes=False)
        t_fin = (time.perf_counter() - t0) / len(specs) / 2
        return {"sec_per_branch_total_derivative_fd_impl": t_der, "sec_per_branch_finite_resolve_eig": t_fin}
    raise ValueError(t)


def aggregate():
    recs = [r for r in I.Store(PHASE).all() if r.get("ok")]
    out = {"label": "POST-HOC exploratory (H31)"}
    part = [r for r in recs if r["task"]["kind"] == "part"]
    grp = {}
    for r in part:
        for k, v in r["groups"].items():
            grp.setdefault(k, []).append(v)
    out["fast_participation_mean_top"] = dict(sorted(((k, float(np.sum(v) / len(part))) for k, v in grp.items()), key=lambda kv: -kv[1])[:10])
    out["fast_lambda_quantiles"] = np.percentile([r["lam"][0] for r in part], [5, 50, 95]).round(1).tolist()
    out["fast_imag_max"] = float(max(abs(r["lam"][1]) for r in part))
    out["fast_sigma_min_gz_quantiles"] = np.percentile([r["sigma_min_gz"] for r in part], [5, 50, 95]).tolist()
    sw = [r for r in recs if r["task"]["kind"] == "sweep"]
    out["sweeps"] = []
    for r in sw:
        rows = [x for x in r["rows"] if "alpha" in x]
        a = np.array([x["alpha"] for x in rows])
        s = np.array([x["sigma_min_gz"] for x in rows])
        rho = np.array([x["rho"] for x in rows])
        cross = rho[np.argmax(a > 0)] if (a > 0).any() else None
        out["sweeps"].append({"pid": r["task"]["pid"], "i": r["task"]["i"], "rho_first_unstable": cross, "alpha_end": float(a[-1]),
                              "sigma_min_gz_min_over_sweep": float(s.min()), "sigma_min_gz_at_crossing": float(s[np.argmax(a > 0)]) if cross is not None else None,
                              "max_alpha_step_jump": float(np.max(np.diff(a))), "n_errors": len(r["rows"]) - len(rows)})
    co = [r for r in recs if r["task"]["kind"] == "cont"]
    out["continuation_n"] = len(co)
    out["continuation_match_frac"] = float(np.mean([r["continuation_matches_crit"] for r in co])) if co else np.nan
    ex = co[:6]
    out["continuation_examples"] = [{"pid": r["task"]["pid"], "S": r["task"]["S"], "i": r["task"]["i"], "match": r["continuation_matches_crit"]} for r in ex]
    tm = [r for r in recs if r["task"]["kind"] == "timing"]
    if tm:
        out["timing"] = {k: v for k, v in tm[0].items() if k.startswith("sec")}
    HI.write_json("H31_explore.json", out)
    return out


if __name__ == "__main__":
    import sys

    if "agg" in sys.argv:
        print(json.dumps(aggregate(), indent=1, default=str))
    else:
        ts = tasks()
        I.run_tasks(PHASE, "H31_explore", "run_task", ts, 6)
        print(json.dumps(aggregate(), indent=1, default=str))
