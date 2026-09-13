# ruff: noqa: E501
"""R10/R11 Model A branch reinforcement on IEEE-68 (target V68 fully converted; REAL).

Per task (policy or draw, branch chunk): SPR-68 IFT derivatives (total, frozen; critical transverse eigenvalue and,
post hoc, the rightmost EM-band mode) and the independent finite re-equilibrated truth at gamma = 1.10/1.25/1.50;
R11 small central steps (gamma = 1 +- 1e-4) for the preregistered branches at the discovery policies;
D_conv: the same derivatives for the all-SG base.

    python R10_branches.py main    -> raw/R10A (V68), raw/R10C (base, D_conv)
    python R10_branches.py draws   -> raw/R15L (ranking subset of the A draws, reference policy)
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import sys
import time
from multiprocessing import Pool

import _sens68 as S

WORKERS = 10
CHUNKS = 4


def design():
    return json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())


def _chunks(lst, n):
    k = (len(lst) + n - 1) // n
    return [lst[i: i + k] for i in range(0, len(lst), k)]


def tasks_main():
    d = design()
    out = []
    r11 = set(d["r11_branches"])
    for p in d["policies_A"]:
        for ch in _chunks(d["eligible_branches"], CHUNKS):
            small = [e for e in ch if e in r11] if p["split"] == "discovery" else []
            out.append({"phase": "R10A", "pid": p["id"], "g": p["g"], "k": p["k"], "target": "V68", "branches": ch, "small": small})
            out.append({"phase": "R10C", "pid": p["id"], "g": p["g"], "k": p["k"], "target": "BASE", "branches": ch, "small": []})
    return out


def tasks_draws():
    import R06_census as C6
    d = design()
    p = C6.reference_policy()
    out = []
    for did in d["draws_ranking_subset_A"]:
        for ch in _chunks(d["eligible_branches"], CHUNKS):
            out.append({"phase": "R15L", "pid": p["id"], "g": p["g"], "k": p["k"], "target": "V68", "branches": ch, "small": [], "draw": did})
    return out


_DRAWS: dict = {}


def _init():
    R.env_threads()
    for dr in design()["draws_A"]:
        _DRAWS[dr["id"]] = {"machine": dr["machine"], "converter": dr["converter"]}


def run_task(task):
    t0 = time.time()
    try:
        members = R.V68 if task["target"] == "V68" else ()
        draw = _DRAWS.get(task.get("draw")) if task.get("draw") else None
        eng = S.Engine68(members, g=task["g"], k=task["k"], draw=draw)
        der = eng.derivatives(task["branches"])
        lam = complex(*der["lam"])
        lam_em = complex(*der["lam_em"]) if der.get("lam_em") else None
        fin = {}
        if task["target"] == "V68":
            for e in task["branches"]:
                gs = [1.10, 1.25, 1.50] + ([1.0001, 0.9999] if e in task["small"] else [])
                fin[str(e)] = {str(gm): S.finite_tracked(members, g=task["g"], k=task["k"], variant="REAL", e=e, gamma=gm, lam_ref=lam, lam_em=lam_em, draw=draw) for gm in gs}
        return task, {"der": der, "finite": fin, "wall_s": time.time() - t0}, True
    except Exception as ex:  # noqa: BLE001
        return task, {"error": repr(ex)[:400]}, False


def run(tasks):
    by_phase = {}
    for t in tasks:
        by_phase.setdefault(t["phase"], []).append(t)
    stores = {ph: R.Store(ph) for ph in by_phase}
    todo = [t for t in tasks if not stores[t["phase"]].done(t)]
    R.log("R10", f"{len(tasks)} tasks, {len(todo)} to run, workers={WORKERS}")
    t0 = time.time()
    with Pool(WORKERS, initializer=_init) as pool:
        for i, (task, rec, ok) in enumerate(pool.imap_unordered(run_task, todo, chunksize=1), 1):
            stores[task["phase"]].put(task, rec, ok)
            R.log("R10", f"{i}/{len(todo)} {task['phase']} {task['pid']} ok={ok} ({time.time() - t0:.0f}s)")
    R.log("R10", f"COMPLETE in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    run(tasks_main() if mode == "main" else tasks_draws())
