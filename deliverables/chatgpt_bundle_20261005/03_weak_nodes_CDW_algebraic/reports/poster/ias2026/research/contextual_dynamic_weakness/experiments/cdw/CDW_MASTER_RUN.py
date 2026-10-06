# ruff: noqa: E501
"""CDW autonomous orchestrator.

Runs the preregistered phases (docs/CDW_PREREG_V1.md) in dependency order with per-task
checkpoints under raw/<phase>/ (resume-safe: completed tasks are never recomputed),
per-phase logs under logs/, atomic status in results/CDW_RUN_STATUS.json, and a run
manifest. A phase whose tasks error is FAILED and blocks only its dependents; a
scientific gate failure is recorded (in the phase gate JSON) and blocks nothing.

Usage: python CDW_MASTER_RUN.py [--only E01,E34] [--max-workers 8]
Environment: BLAS pinned to one thread (see _infra). Never touches TX4 files.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import _infra as I  # noqa: E402

# (name, module, tasks fn, run fn, aggregate fn, dependencies, pre-hook)
PHASES = [
    ("E01", "E01_census", "tasks", "run_task", "aggregate", [], None),
    ("E34", "E34_sens", "tasks", "run_task", "aggregate", [], None),
    ("E09", "E09_design", "tasks", "run_task", "aggregate", [], None),
    ("E02", "E02_robust", "tasks", "run_task", "aggregate", [], None),
    ("E06", "E06_corridors", "tasks", "run_task", "aggregate", ["E34"], None),
    ("E08a", "E08_exchange", "tasks_a", "run_task_a", None, [], None),
    ("E08b", "E08_exchange", "tasks_b", "run_task_b", "aggregate", ["E08a"], None),
    ("E12", "E11_spectral", "tasks12", "run_task12", None, [], None),
    ("E12q3", "E11_spectral", "tasks12q3", "run_task12q3", "aggregate12", ["E01", "E12"], None),
    ("E13", "E13_reduction", "tasks13", "run_task13", "aggregate13", ["E01"], "build_pod_basis"),
    ("E14", "E13_reduction", "tasks14", "run_task14", "aggregate14", [], None),
    ("GC", "E13_reduction", None, None, "gold_c", ["E13", "E14"], None),
    ("E16", "E16_modal_energy", "tasks", "run_task", "aggregate", [], None),
    ("E17", "E17_africano_gate", None, None, "run", [], None),
    ("E07", "E07_topology", "tasks", "run_task", "aggregate", [], None),
    ("E23c", "E23_crossmodel", "tasks_c", "run_task_c", None, [], None),
    ("E23l", "E23_crossmodel", "tasks_l", "run_task_l", "aggregate", ["E01", "E34", "E06", "E23c"], None),
    ("E11", "E11_spectral", None, None, "aggregate11", ["E01", "E34"], None),
    ("E05", "E05_static_dynamic", None, None, "aggregate", ["E01", "E11", "E34"], None),
]


class Log:
    def __init__(self, phase):
        self.f = (I.LOGS / f"{phase}.log").open("a", encoding="utf-8")

    def __call__(self, msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, flush=True)
        self.f.write(line + "\n")
        self.f.flush()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--max-workers", type=int, default=12)
    args = ap.parse_args(argv)
    only = {x for x in args.only.split(",") if x}
    st0 = I.load_status()
    st0["environment"] = I.environment()
    st0["started"] = st0.get("started") or time.strftime("%Y-%m-%dT%H:%M:%S")
    I.atomic_write_json(I.STATUS, st0)
    manifest = []
    for name, module, tfn, rfn, afn, deps, pre in PHASES:
        if only and name not in only:
            continue
        log = Log(name)
        st = I.load_status()["phases"]
        if st.get(name, {}).get("state") == "COMPLETE" and not only:
            log(f"[{name}] already COMPLETE, skipping")
            continue
        bad = [d for d in deps if I.load_status()["phases"].get(d, {}).get("state") in ("FAILED", "BLOCKED")]
        if bad:
            I.set_status(name, state="BLOCKED", reason=f"dependency failed: {bad}")
            log(f"[{name}] BLOCKED by {bad}")
            continue
        t0 = time.time()
        res = I.resources(args.max_workers)
        I.set_status(name, state="RUNNING", started=time.strftime("%Y-%m-%dT%H:%M:%S"), resources=res)
        try:
            mod = importlib.import_module(module)
            if pre:
                log(f"[{name}] pre-hook {pre}: {getattr(mod, pre)()}")
            info = {}
            if tfn:
                tasks = getattr(mod, tfn)()
                info = I.run_tasks(name, module, rfn, tasks, res["workers"], log=log)
                if info["n_tasks"] and _n_errors(name) > 0.02 * info["n_tasks"]:
                    raise RuntimeError(f"{_n_errors(name)} task errors (>2%)")
            gate = getattr(mod, afn)() if afn else None
            if hasattr(gate, "to_dict"):
                gate = gate.to_dict("records")
            I.set_status(name, state="COMPLETE", wall_s=round(time.time() - t0, 1), tasks=info, n_errors=_n_errors(name),
                         gate=_compact(gate))
            log(f"[{name}] COMPLETE in {time.time() - t0:.0f}s; gate={json.dumps(_compact(gate), default=str)[:600]}")
        except Exception as e:  # noqa: BLE001
            I.set_status(name, state="FAILED", wall_s=round(time.time() - t0, 1), error=f"{type(e).__name__}: {e}",
                         trace=traceback.format_exc()[-3000:])
            log(f"[{name}] FAILED: {e}\n{traceback.format_exc()}")
        manifest.append({"phase": name, "module": module, "state": I.load_status()["phases"][name]["state"],
                         "wall_s": round(time.time() - t0, 1)})
    st = I.load_status()
    st["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    I.atomic_write_json(I.STATUS, st)
    return 0


def _n_errors(phase):
    return sum(1 for p in (I.RAW / phase).glob("*.json") if '"ok": false' in p.read_text(encoding="utf-8")[:40])


def _compact(g):
    try:
        s = json.dumps(g, default=str)
    except Exception:  # noqa: BLE001
        return str(g)[:2000]
    return json.loads(s) if len(s) < 20000 else s[:2000]


if __name__ == "__main__":
    raise SystemExit(main())
