# ruff: noqa: E402
"""Run PD2 (24 design tasks; per-task JSON checkpoints, resumable)."""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

import sys
import time
from multiprocessing import Pool

import pd2_census_design as PD2


def main() -> int:
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    tasks = PD2.tasks()
    print(f"PD2: {len(tasks)} tasks on {workers} workers", flush=True)
    t0 = time.time()
    with Pool(processes=workers) as pool:
        for rec in pool.imap_unordered(PD2.run_and_store, tasks, chunksize=1):
            t = rec.get("task", {})
            print(f"  {PD2.task_name(t) if t else '?'}: ok={rec.get('ok')} stop={rec.get('stop')} "
                  f"phi {rec.get('phi_nominal')} -> {rec.get('phi_final')}  {time.time() - t0:7.0f} s", flush=True)
    print(f"PD2 finished in {time.time() - t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
