# ruff: noqa: E402
"""Run PD1 over the 240 unstable E35 samples (checkpointed JSONL, resumable)."""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

import json
import sys
import time
from multiprocessing import Pool

import _pd_common as P
import pandas as pd

import pd1_adaptive_retune as PD1


def samples() -> list[dict]:
    t = pd.read_parquet(P.E35_TABLE)
    u = t[(t.status == "ACCEPTED") & (t.portfolio_unstable == True)]  # noqa: E712
    cols = ["sample", "load", "reactive_load", "availability", "alpha_IA_exact", "rc_success", "sc_success"]
    return u[cols].to_dict("records")


def main() -> int:
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    out = P.out_dir("PD1") / "PD1_records.jsonl"
    done = set()
    if out.exists():
        done = {json.loads(line)["sample"] for line in out.read_text(encoding="utf-8").splitlines() if line.strip()}
    todo = [r for r in samples() if int(r["sample"]) not in done]
    print(f"PD1: {len(todo)} samples to run on {workers} workers ({len(done)} done)", flush=True)
    t0 = time.time()
    with Pool(processes=workers, initializer=PD1.initialise) as pool, out.open("a", encoding="utf-8") as fh:
        for k, rec in enumerate(pool.imap_unordered(PD1.run_sample, todo, chunksize=1), 1):
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            if k % 10 == 0:
                print(f"  {k}/{len(todo)}  {time.time() - t0:7.0f} s", flush=True)
    print(f"PD1 finished in {time.time() - t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
