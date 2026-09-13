# ruff: noqa: E402, E501
"""PD2 - census-scale plan-level design vs single-boundary tuning (re-preregistered E9).

Preregistration: docs/20260913_PLANNING_DESIGN_PREREG.md.

Unchanged from E9 (CDW prereg b8ae3082): targets, policies, parameter families,
bounds, ranges, the linearized QP step (cvxpy CLARABEL/SCS/OSQP with the Gordan
least-violation fallback), EPS, the stopping rules and the GOLD-D criterion.
Changed, fixed before any PD2 run:
  - iteration budget 6 -> 30, with a cap of 40 lattice evaluations and 6 h per task;
  - step acceptance: a step is accepted only if the merit decreases (merit = alpha_T
    for "single", Phi = max over the lattice for "plan"); otherwise the trust radius
    is halved (floor TRUST/16) and the QP is re-solved with the same gradients;
  - plan active set: every subset with alpha >= -EPS, ranked by alpha, cap 60
    (E9: subsets within 0.05 of Phi, cap 25); subsets with non-finite alpha are skipped.
"""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

import json
import sys
import time

import _pd_common as P
import numpy as np

sys.path.insert(0, str(P.CDW_EXPERIMENTS))
import _analysis as AN
import _cdw as C
import E09_design as E9

TARGETS = {"T1": ("H4", "D01"), "T2": ("V9", "D01"), "T3": ("H4", "D11"), "T4": ("V9", "D03")}
SETS = {"H4": C.V4, "V9": C.V9}
FAMILIES = E9.FAMILIES
METHODS = E9.METHODS
EPS = E9.EPS
BOUNDS, RANGE, TRUST0 = E9.BOUNDS, E9.RANGE, E9.TRUST
ITERS = 30
MAX_LATTICE_EVALS = 40
MAX_WALL_S = 6 * 3600.0
MAX_ACTIVE = 60
TRUST_FLOOR = TRUST0 / 16.0
TASK_DIR = P.RESULTS / "PD2" / "tasks"


def tasks():
    # long tasks first (plan at census scale)
    order = [(t, f, m) for t in ("T2", "T4") for f in FAMILIES for m in ("plan", "single")]
    order += [(t, f, m) for t in ("T1", "T3") for f in FAMILIES for m in METHODS]
    return [{"target": t, "family": f, "method": m} for t, f, m in order]


def task_name(task) -> str:
    return f"{task['target']}_{task['family']}_{task['method']}"


def qp_step(active, alphas, grads, vals, plist, trust):
    """E9.qp_step with an explicit trust radius (everything else verbatim)."""

    import cvxpy as cp

    n = len(plist)
    a = np.array([vals[p] for p in plist])
    lo = np.array([BOUNDS[k][0] for k, _ in plist]) - a
    hi = np.array([BOUNDS[k][1] for k, _ in plist]) - a
    rng = np.array([RANGE[k] for k, _ in plist])
    tr = trust * rng
    lo, hi = np.maximum(lo, -tr), np.minimum(hi, tr)
    d = cp.Variable(n)
    G = np.array([grads[s] for s in active])
    al = np.array([alphas[s] for s in active])
    prob = cp.Problem(cp.Minimize(0.5 * cp.sum_squares(cp.multiply(d, 1 / rng))), [d >= lo, d <= hi, al + G @ d <= -EPS])
    status = None
    for solver in (cp.CLARABEL, cp.SCS, cp.OSQP):
        try:
            prob.solve(solver=solver)
            status = prob.status
            if status in ("optimal", "optimal_inaccurate"):
                break
        except Exception:  # noqa: BLE001
            continue
    if status in ("optimal", "optimal_inaccurate"):
        return np.asarray(d.value).ravel(), "feasible", None
    gv, lam = AN.gordan(G) if len(active) > 1 else (float(np.sum(G**2)), np.ones(1))
    t = cp.Variable()
    d2 = cp.Variable(n)
    prob2 = cp.Problem(cp.Minimize(t + 1e-3 * cp.sum_squares(cp.multiply(d2, 1 / rng))), [d2 >= lo, d2 <= hi, al + G @ d2 <= t])
    st2 = None
    for solver in (cp.CLARABEL, cp.SCS, cp.OSQP):
        try:
            prob2.solve(solver=solver)
            st2 = prob2.status
            if st2 in ("optimal", "optimal_inaccurate") and d2.value is not None:
                break
        except Exception:  # noqa: BLE001
            continue
    if st2 in ("optimal", "optimal_inaccurate") and d2.value is not None:
        return np.asarray(d2.value).ravel(), "infeasible", {"gordan_value": gv}
    return np.zeros(n), "solver_failed", {"gordan_value": gv}


def _finite(x) -> float:
    return float(x) if np.isfinite(x) else np.inf


def run_task(task: dict) -> dict:
    TASK_DIR.mkdir(parents=True, exist_ok=True)
    tname, pid = TARGETS[task["target"]]
    T = tuple(SETS[tname])
    theta = C.policy(pid)
    plist = E9.params(T, task["family"])
    vals = {p: (theta[0] if p[0] == "g" else 1.0) for p in plist}
    t0 = time.time()
    n_lat = 0

    def lattice(v):
        nonlocal n_lat
        n_lat += 1
        return E9.evaluate_all(T, theta, v)

    def merit_single(v):
        r = C.eval_portfolio(T, theta, modes=False, **E9.S.design_kwargs(E9.design_of(v)))
        return _finite(r["alpha"])

    ev = lattice(vals)
    phi0 = max(_finite(a) for a, _ in ev.values())
    rec = {"task": task, "phi_nominal": phi0, "alpha_T_nominal": ev[T][0],
           "n_unstable_nominal": sum(st == "UNSTABLE" for _, st in ev.values())}
    if phi0 < 0:
        rec["outcome"] = "ALREADY_SAFE"
        return rec
    hist, active, trust = [], set(), TRUST0
    alphas = {s: a for s, (a, _) in ev.items()}
    aT = _finite(alphas[T])
    phi = phi0
    stop = "ITERS"
    for it in range(ITERS):
        done = (aT <= -EPS / 2) if task["method"] == "single" else (phi <= -EPS / 2)
        h = {"iter": it, "phi": phi, "alpha_T": aT, "trust": trust, "lattice_evals": n_lat, "wall_s": time.time() - t0}
        if done:
            hist.append(h)
            stop = "CONVERGED"
            break
        if n_lat >= MAX_LATTICE_EVALS or time.time() - t0 > MAX_WALL_S:
            hist.append(h)
            stop = "BUDGET_EXHAUSTED"
            break
        if task["method"] == "single":
            active = {T}
        else:
            cand = sorted((s for s, a in alphas.items() if np.isfinite(a) and a >= -EPS), key=lambda s: -alphas[s])
            keep = {s for s in active if np.isfinite(alphas.get(s, np.nan)) and alphas[s] >= -EPS}
            active = set(sorted(set(cand[:MAX_ACTIVE]) | keep, key=lambda s: -alphas[s])[:MAX_ACTIVE])
        grads = {s: E9.gradient(s, theta, vals, plist) for s in active}
        merit_old = aT if task["method"] == "single" else phi
        accepted = False
        tries = []
        while trust >= TRUST_FLOOR - 1e-15:
            d, qst, conflict = qp_step(sorted(active), alphas, grads, vals, plist, trust)
            trial = {p: float(min(max(vals[p] + d[j], BOUNDS[p[0]][0]), BOUNDS[p[0]][1])) for j, p in enumerate(plist)}
            if task["method"] == "single":
                m_new = merit_single(trial)
                tries.append({"trust": trust, "qp": qst, "merit": m_new})
                if m_new < merit_old:
                    vals, aT, accepted = trial, m_new, True
                    break
            else:
                ev_t = lattice(trial)
                al_t = {s: a for s, (a, _) in ev_t.items()}
                m_new = max(_finite(a) for a in al_t.values())
                tries.append({"trust": trust, "qp": qst, "merit": m_new, "gordan": (conflict or {}).get("gordan_value")})
                if m_new < merit_old:
                    vals, ev, alphas, phi, accepted = trial, ev_t, al_t, m_new, True
                    aT = _finite(alphas[T])
                    break
                if n_lat >= MAX_LATTICE_EVALS:
                    break
            trust /= 2.0
        h.update(n_active=len(active), tries=tries, accepted=accepted)
        hist.append(h)
        if not accepted:
            stop = "NO_DESCENT" if trust < TRUST_FLOOR else "BUDGET_EXHAUSTED"
            break
        trust = min(TRUST0, trust * 2.0)  # expand after an accepted step
        P.write_json(TASK_DIR / f"{task_name(task)}.partial.json", {**rec, "history": hist})
    ev = lattice(vals) if task["method"] == "single" else ev
    alphas = {s: a for s, (a, _) in ev.items()}
    statuses = {s: st for s, (_, st) in ev.items()}
    rec.update(
        outcome="DESIGNED", stop=stop, history=hist,
        phi_final=max(_finite(a) for a in alphas.values()), alpha_T_final=alphas[T], status_T_final=statuses[T],
        unstable_subsets_final=[C.label(s) for s, st in statuses.items() if st == "UNSTABLE"],
        n_unresolved_final=sum(st == "BOUNDARY_OR_UNRESOLVED" for st in statuses.values()),
        n_infeasible_final=sum(st == "INFEASIBLE" for st in statuses.values()),
        final_values={f"{k}{i}": v for (k, i), v in vals.items() if abs(v - (theta[0] if k == "g" else 1.0)) > 1e-9},
        lattice_evals=n_lat, wall_s=time.time() - t0,
    )
    return rec


def run_and_store(task: dict) -> dict:
    TASK_DIR.mkdir(parents=True, exist_ok=True)
    path = TASK_DIR / f"{task_name(task)}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    try:
        rec = run_task(task)
        rec["ok"] = True
    except Exception as e:  # noqa: BLE001 - recorded, never silently dropped
        rec = {"task": task, "ok": False, "error": repr(e)[:400]}
    P.write_json(path, rec)
    return rec


def gold_d(records: list[dict]) -> dict:
    """E9 GOLD-D criterion, verbatim."""

    by = {(r["task"]["target"], r["task"]["family"], r["task"]["method"]): r for r in records if r.get("ok") and r.get("outcome") == "DESIGNED"}
    cases, gold = [], False
    for (tg, fam, m), s in by.items():
        if m != "single" or (tg, fam, "plan") not in by:
            continue
        p = by[(tg, fam, "plan")]
        single_final_only = (s["status_T_final"] == "STABLE") and (s["phi_final"] >= 0)
        plan_safe = (p["phi_final"] < 0) and (p["n_unresolved_final"] == 0) and (not p["unstable_subsets_final"])
        cases.append({"target": tg, "family": fam, "single_T_stable": s["status_T_final"] == "STABLE", "single_phi": s["phi_final"],
                      "single_leaves_unsafe_subset": bool(single_final_only), "plan_phi": p["phi_final"], "plan_safe": bool(plan_safe),
                      "single_stop": s.get("stop"), "plan_stop": p.get("stop")})
        gold |= bool(single_final_only and plan_safe)
    return {"GOLD_D_pass": gold, "cases": cases}
