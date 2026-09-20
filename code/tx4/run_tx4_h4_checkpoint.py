"""Checkpointed exact H4 pass for the declared TX4 condition design."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tx4_robustness as exact  # noqa: E402

FIELDS = [*exact.FIELDS, "evaluator_tier"]
H4 = (30, 33, 35, 37)


def work(item):
    campaign, condition_id, seed, params = item
    row = exact._one_portfolio(params, H4, campaign, condition_id, seed)
    row["evaluator_tier"] = "EXACT_H4_DAE"
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=max(1, min(18, (os.cpu_count() or 2) - 2)))
    args = parser.parse_args()
    conditions = exact.make_conditions()
    calibration = set()
    for item in conditions:
        cid = item[1]
        if cid == "nominal":
            calibration.add(cid)
        elif cid.startswith("1d_"):
            idx = int(cid.rsplit("_", 1)[-1])
            if idx % 10 == 0 or idx == 200:
                calibration.add(cid)
        elif cid.startswith("qmc_") and int(cid.split("_")[1]) < 64:
            calibration.add(cid)
        elif cid.startswith("mc_") and int(cid.split("_")[1]) < 64:
            calibration.add(cid)
        elif cid.startswith("2d_"):
            parts = cid.rsplit("_", 2)
            if int(parts[-2]) % 10 == 0 and int(parts[-1]) % 10 == 0:
                calibration.add(cid)
    remaining = [item for item in conditions if item[1] not in calibration]
    batch = remaining[args.start : args.start + args.count]
    out = ROOT / "results" / "TX4_ROBUSTNESS_H4_CHECKPOINTS.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    exists = out.exists()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(work, item) for item in batch]
        rows = [future.result() for future in as_completed(futures)]
    with out.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
    print(f"TX4_H4_CHECKPOINT start={args.start} count={len(rows)} total_remaining={len(remaining)}")


if __name__ == "__main__":
    main()
