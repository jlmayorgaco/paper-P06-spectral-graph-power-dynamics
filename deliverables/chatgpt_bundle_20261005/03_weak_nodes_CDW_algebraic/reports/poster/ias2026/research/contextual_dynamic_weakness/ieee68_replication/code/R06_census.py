# ruff: noqa: E501
"""R6 exhaustive portfolio census, Model A (TX4 GFL): REAL x 16 policies, NOGOV/SP33 x 6 ablation policies,
and (phase R15A) REAL x 40 envelope draws at the reference policy. Resume-safe (one JSON per task).

    python R06_census.py main      -> raw/R06A/*.json
    python R06_census.py draws     -> raw/R15A/*.json   (reference policy fixed by R15 rule)
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import sys
import time
from multiprocessing import Pool

WORKERS = 16


def design():
    return json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())


def tasks_main():
    d = design()
    pol = {p["id"]: p for p in d["policies_A"]}
    out = []
    for p in d["policies_A"]:
        for S in R.subsets():
            out.append({"phase": "R06A", "variant": "REAL", "pid": p["id"], "g": p["g"], "k": p["k"], "S": R.label(S)})
    for pid in d["ablation_policies"]:
        p = pol[pid]
        for v in ("NOGOV", "SP33"):
            for S in R.subsets():
                out.append({"phase": "R06A", "variant": v, "pid": pid, "g": p["g"], "k": p["k"], "S": R.label(S)})
    return out


def reference_policy():
    """R15 rule: base-stable policy nearest (u=0.5, k=1.25) in normalized coordinates; ties -> lower index."""
    d = design()
    st = R.Store("R06A")
    base = {}
    for r in st.all():
        t = r["task"]
        if t["variant"] == "REAL" and t["S"] == "BASE" and r.get("ok"):
            base[t["pid"]] = r["record"]["status"]
    cands = []
    for i, p in enumerate(d["policies_A"]):
        if base.get(p["id"]) == "STABLE":
            du = (p["u"] - 0.5) / 1.0
            dk = (p["k"] - 1.25) / 1.5
            cands.append(((du * du + dk * dk) ** 0.5, i, p))
    assert cands, "no base-stable policy"
    return sorted(cands, key=lambda c: (c[0], c[1]))[0][2]


def tasks_draws():
    d = design()
    p = reference_policy()
    out = []
    for dr in d["draws_A"]:
        for S in R.subsets():
            out.append({"phase": "R15A", "variant": "REAL", "pid": p["id"], "g": p["g"], "k": p["k"], "S": R.label(S), "draw": dr["id"]})
    return out


_DRAWS: dict = {}


def _init():
    R.env_threads()
    for dr in design()["draws_A"]:
        _DRAWS[dr["id"]] = {"machine": dr["machine"], "converter": dr["converter"]}


def run_task(task):
    members = () if task["S"] == "BASE" else tuple(int(b) for b in task["S"].split("+"))
    draw = _DRAWS.get(task.get("draw")) if task.get("draw") else None
    try:
        rec = R.eval_portfolio(members, variant=task["variant"], g=task["g"], k=task["k"], draw=draw)
        return task, rec, True
    except Exception as e:  # noqa: BLE001
        return task, {"error": repr(e)[:300]}, False


def run(phase_tasks, phase):
    st = R.Store(phase)
    todo = [t for t in phase_tasks if not st.done(t)]
    R.log(phase, f"{len(phase_tasks)} tasks, {len(phase_tasks) - len(todo)} done, {len(todo)} to run, workers={WORKERS}")
    t0 = time.time()
    n_ok = n_err = 0
    with Pool(WORKERS, initializer=_init) as pool:
        for i, (task, rec, ok) in enumerate(pool.imap_unordered(run_task, todo, chunksize=2), 1):
            st.put(task, rec, ok)
            n_ok += ok
            n_err += not ok
            if i % 100 == 0:
                R.log(phase, f"{i}/{len(todo)} ({time.time() - t0:.0f}s)")
    R.log(phase, f"COMPLETE ok={n_ok} err={n_err} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    if mode == "main":
        run(tasks_main(), "R06A")
    elif mode == "draws":
        run(tasks_draws(), "R15A")
