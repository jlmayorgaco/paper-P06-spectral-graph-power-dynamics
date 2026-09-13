# ruff: noqa: E501
"""CDW hardening orchestrator (compute phases). Prereg: docs/CDW_HARDENING_PREREG_V1.md.

Per-task checkpoints under raw/H_<phase>/ (resume-safe: completed tasks are never recomputed),
per-phase logs under logs/hardening/, status in results/hardening/CDWH_RUN_STATUS.json.
Analyses (H03-H05, H07, H11, H13, H15, H16, H19) run separately from the stored raw data.

Usage: python CDWH_MASTER_RUN.py [--only H01,H06] [--max-workers 20]
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE.parent / "cdw")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import _infra as I  # noqa: E402  (BLAS pinned before numpy)
import _hinfra as HI  # noqa: E402

# (name, store, module, tasks fn, run fn, dependencies)
PHASES = [
    ("H01", "H_H01", "H01_census", "tasks", "run_task", []),
    ("H01e", "H_H01e", "H01_census", "tasks_e", "run_task_e", []),
    ("H06", "H_H06", "H06_links", "tasks", "run_task", []),
    ("H08n", "H_H08n", "H06_links", "tasks_node", "run_task", []),
    ("H09", "H_H09", "H09_corridors", "tasks9", "run_task9", []),
    ("H12", "H_H12", "H12_mixing", "tasks", "run_task", []),
    ("H12b", "H_H12b", "H12_mixing", "tasks_b", "run_task_b", ["H12"]),
    ("H14", "H_H14", "H14_topology", "tasks", "run_task", []),
    ("H18c", "H_H18c", "H14_topology", "tasks_c", "run_task_c", []),
    ("H10", "H_H10", "H09_corridors", "tasks10", "run_task10", []),
]


class Log:
    def __init__(self, phase):
        self.f = (HI.LOGS / f"{phase}.log").open("a", encoding="utf-8")

    def __call__(self, msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, flush=True)
        self.f.write(line + "\n")
        self.f.flush()


def n_errors(store):
    return sum(1 for p in (I.RAW / store).glob("*.json") if '"ok": false' in p.read_text(encoding="utf-8")[:40])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--max-workers", type=int, default=20)
    args = ap.parse_args(argv)
    only = [x for x in args.only.split(",") if x]
    st = HI.load_status()
    st["environment"] = I.environment()
    st["prereg_commit"] = HI.PREREG_COMMIT
    I.atomic_write_json(HI.STATUS, st)
    for name, store, module, tfn, rfn, deps in PHASES:
        if only and name not in only:
            continue
        log = Log(name)
        bad = [d for d in deps if HI.load_status()["phases"].get(d, {}).get("state") != "COMPLETE"]
        if bad:
            HI.set_status(name, state="BLOCKED", reason=f"dependency not complete: {bad}")
            log(f"[{name}] BLOCKED by {bad}")
            continue
        t0 = time.time()
        res = I.resources(args.max_workers)
        HI.set_status(name, state="RUNNING", started=time.strftime("%Y-%m-%dT%H:%M:%S"), resources=res)
        try:
            mod = importlib.import_module(module)
            tasks = getattr(mod, tfn)()
            info = I.run_tasks(store, module, rfn, tasks, res["workers"], log=log)
            ne = n_errors(store)
            state = "COMPLETE" if ne <= 0.02 * max(info["n_tasks"], 1) else "FAILED"
            HI.set_status(name, state=state, wall_s=round(time.time() - t0, 1), tasks=info, n_errors=ne)
            log(f"[{name}] {state} in {time.time() - t0:.0f}s; errors={ne}")
        except Exception as e:  # noqa: BLE001
            HI.set_status(name, state="FAILED", wall_s=round(time.time() - t0, 1), error=f"{type(e).__name__}: {e}",
                          trace=traceback.format_exc()[-3000:])
            log(f"[{name}] FAILED: {e}\n{traceback.format_exc()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
