"""CDW infrastructure: paths, atomic writes, checkpoints, resource policy, pool, manifest.

Nothing here touches TX4 files. All writes go under contextual_dynamic_weakness/.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# BLAS pinning must happen before numpy is imported anywhere in the process.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]  # contextual_dynamic_weakness/
RESEARCH = PROJECT.parent  # reports/poster/ias2026/research/
REPO = RESEARCH.parents[3]
RESULTS = PROJECT / "results"
RAW = PROJECT / "raw"
LOGS = PROJECT / "logs"
FIGURES = PROJECT / "figures"
DOCS = PROJECT / "docs"
STATUS = RESULTS / "CDW_RUN_STATUS.json"
PREREG_COMMIT = "b8ae3082"

for p in (RESULTS, RAW, LOGS, FIGURES):
    p.mkdir(parents=True, exist_ok=True)


def atomic_write_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def atomic_write_json(path: Path, payload) -> None:
    atomic_write_text(path, json.dumps(payload, indent=1, default=_default, sort_keys=False))


def _default(o):
    import numpy as np

    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def task_key(task: dict) -> str:
    return hashlib.sha1(json.dumps(task, sort_keys=True, default=str).encode()).hexdigest()[:16]


# ------------------------------------------------------------------ resources --
def resources(max_workers: int | None = None, mem_per_worker_gb: float = 0.35) -> dict:
    import psutil

    ncpu = os.cpu_count() or 2
    vm = psutil.virtual_memory()
    free_gb = vm.available / 2**30
    by_cpu = max(1, ncpu - 2)
    by_mem = max(1, int((free_gb - 1.5) / mem_per_worker_gb))
    w = min(by_cpu, by_mem, max_workers or 12)
    return {"cpu": ncpu, "free_gb": round(free_gb, 2), "total_gb": round(vm.total / 2**30, 2), "workers": w}


def git_head() -> str:
    try:
        return subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def environment() -> dict:
    import numpy
    import pandas
    import scipy

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "git_head": git_head(),
        "prereg_commit": PREREG_COMMIT,
        "blas_threads": {v: os.environ.get(v) for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
    }


# --------------------------------------------------------------- checkpoints --
class Store:
    """One JSON file per completed task under raw/<phase>/; resume skips done tasks."""

    def __init__(self, phase: str):
        self.dir = RAW / phase
        self.dir.mkdir(parents=True, exist_ok=True)

    def path(self, key: str) -> Path:
        return self.dir / f"{key}.json"

    def done(self, key: str) -> bool:
        return self.path(key).exists()

    def put(self, key: str, record: dict) -> None:
        atomic_write_json(self.path(key), record)

    def get(self, key: str) -> dict:
        return json.loads(self.path(key).read_text(encoding="utf-8"))

    def all(self) -> list[dict]:
        return [json.loads(p.read_text(encoding="utf-8")) for p in sorted(self.dir.glob("*.json"))]


def _worker_init():
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"


def _safe_call(fn_path: tuple[str, str], task: dict):
    import importlib

    mod = importlib.import_module(fn_path[0])
    fn = getattr(mod, fn_path[1])
    t0 = time.perf_counter()
    try:
        out = fn(task)
        out = {"ok": True, "task": task, **out}
    except Exception as e:  # noqa: BLE001
        out = {"ok": False, "task": task, "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-2000:]}
    out["wall_s"] = round(time.perf_counter() - t0, 3)
    return out


def run_tasks(phase: str, module: str, func: str, tasks: list[dict], workers: int, log=print) -> dict:
    """Run tasks in a process pool with per-task checkpoints. Deterministic per task."""

    store = Store(phase)
    todo = [t for t in tasks if not store.done(task_key(t))]
    log(f"[{phase}] {len(tasks)} tasks, {len(tasks) - len(todo)} already done, {len(todo)} to run, workers={workers}")
    n_ok = n_err = 0
    t0 = time.time()
    if todo:
        if workers <= 1:
            for i, t in enumerate(todo):
                r = _safe_call((module, func), t)
                store.put(task_key(t), r)
                n_ok += r["ok"]
                n_err += not r["ok"]
                if (i + 1) % 50 == 0:
                    log(f"[{phase}] {i + 1}/{len(todo)} ({time.time() - t0:.0f}s)")
        else:
            with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init) as ex:
                futs = {ex.submit(_safe_call, (module, func), t): t for t in todo}
                for i, fut in enumerate(as_completed(futs)):
                    t = futs[fut]
                    try:
                        r = fut.result()
                    except Exception as e:  # noqa: BLE001 (worker crash)
                        r = {"ok": False, "task": t, "error": f"worker: {e}"}
                    store.put(task_key(t), r)
                    n_ok += r["ok"]
                    n_err += not r["ok"]
                    if (i + 1) % 100 == 0:
                        log(f"[{phase}] {i + 1}/{len(todo)} ({time.time() - t0:.0f}s)")
    log(f"[{phase}] finished: ok={n_ok} err={n_err} new, {time.time() - t0:.0f}s")
    return {"n_tasks": len(tasks), "new_ok": n_ok, "new_err": n_err}


# ------------------------------------------------------------------- status --
def load_status() -> dict:
    if STATUS.exists():
        return json.loads(STATUS.read_text(encoding="utf-8"))
    return {"phases": {}, "created": time.strftime("%Y-%m-%dT%H:%M:%S")}


def set_status(phase: str, **fields) -> None:
    st = load_status()
    rec = st["phases"].get(phase, {})
    rec.update(fields)
    rec["updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    st["phases"][phase] = rec
    atomic_write_json(STATUS, st)
