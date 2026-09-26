"""Checkpointed exact all-16 QMC/MC fallback for the TX4 statistics audit."""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tx4_robustness as exact  # noqa: E402


H4_LABEL = "30+33+35+37"
PARAMS = tuple(exact.PARAM_NAMES)
OUTPUT_FIELDS = [*exact.FIELDS, "evaluator_tier", "evaluation_source", "audit_exact"]
CHUNK_SIZE = 250


def row_for_portfolio(item: tuple[str, str, int, dict[str, float], tuple[int, ...]]) -> dict[str, object]:
    campaign, condition_id, seed, params, portfolio = item
    row = exact._one_portfolio(params, portfolio, campaign, condition_id, seed)
    row["evaluator_tier"] = "EXACT_16_AUDIT_DAE"
    row["evaluation_source"] = "EXACT_DAE"
    row["audit_exact"] = True
    return row


def evaluate_condition(item: tuple[str, str, int, dict[str, float]]) -> list[dict[str, object]]:
    campaign, condition_id, seed, params = item
    return [row_for_portfolio((campaign, condition_id, seed, params, portfolio)) for portfolio in exact.PORTFOLIOS]


def conditions_from_master(campaign: str) -> list[tuple[str, str, int, dict[str, float]]]:
    master = pd.read_parquet(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet")
    subset = master[master.campaign == campaign]
    first = subset.drop_duplicates("condition_id", keep="first")
    expected = 4096 if campaign == "QMC" else 5000
    if len(first) != expected:
        raise RuntimeError(f"{campaign}: expected {expected} coordinates, found {len(first)}")
    conditions = []
    for row in first.itertuples(index=False):
        conditions.append(
            (
                str(row.campaign),
                str(row.condition_id),
                int(row.seed),
                {name: float(getattr(row, name)) for name in PARAMS},
            )
        )
    return conditions


def atomic_json(path: Path, value: dict[str, object]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def write_shard(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def run_campaign(campaign: str, workers: int) -> None:
    conditions = conditions_from_master(campaign)
    shard_dir = RESULTS / "TX4_AUDIT_EXACT_SHARDS" / campaign
    shard_dir.mkdir(parents=True, exist_ok=True)
    progress_path = RESULTS / "TX4_AUDIT_EXACT_PROGRESS.json"
    start = time.perf_counter()
    n_chunks = (len(conditions) + CHUNK_SIZE - 1) // CHUNK_SIZE
    completed_chunks = 0
    for chunk_index in range(n_chunks):
        chunk = conditions[chunk_index * CHUNK_SIZE : (chunk_index + 1) * CHUNK_SIZE]
        shard_path = shard_dir / f"shard_{chunk_index:04d}.csv"
        if shard_path.exists():
            completed_chunks += 1
            continue
        rows_by_condition: dict[str, list[dict[str, object]]] = {}
        failures: list[tuple[str, str, int, dict[str, float]]] = []
        try:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(evaluate_condition, item): item for item in chunk}
                for done_index, future in enumerate(as_completed(futures), 1):
                    item = futures[future]
                    try:
                        group = future.result()
                        rows_by_condition[item[1]] = group
                    except Exception as exc:  # retry only the failed condition synchronously
                        print(f"TX4_AUDIT_WORKER_RETRY campaign={campaign} condition={item[1]} error={type(exc).__name__}:{exc}", flush=True)
                        failures.append(item)
                    if done_index % 25 == 0 or done_index == len(chunk):
                        print(
                            f"TX4_AUDIT_PROGRESS campaign={campaign} chunk={chunk_index + 1}/{n_chunks} "
                            f"conditions={done_index}/{len(chunk)}",
                            flush=True,
                        )
        except Exception as exc:
            print(f"TX4_AUDIT_POOL_RETRY campaign={campaign} chunk={chunk_index} error={type(exc).__name__}:{exc}", flush=True)
        failures.extend(item for item in chunk if item[1] not in rows_by_condition and item not in failures)
        for item in failures:
            rows_by_condition[item[1]] = evaluate_condition(item)
        if set(rows_by_condition) != {item[1] for item in chunk}:
            missing = sorted({item[1] for item in chunk} - set(rows_by_condition))
            raise RuntimeError(f"Missing exact conditions after retry: {missing[:5]}")
        shard_rows = [row for item in chunk for row in rows_by_condition[item[1]]]
        if len(shard_rows) != len(chunk) * 16:
            raise RuntimeError(f"{campaign} chunk {chunk_index}: expected {len(chunk) * 16} rows, got {len(shard_rows)}")
        write_shard(shard_path, shard_rows)
        completed_chunks += 1
        atomic_json(
            progress_path,
            {
                "status": "RUNNING",
                "campaign": campaign,
                "completed_chunks": completed_chunks,
                "total_chunks": n_chunks,
                "completed_conditions": min(completed_chunks * CHUNK_SIZE, len(conditions)),
                "total_conditions": len(conditions),
                "chunk_size": CHUNK_SIZE,
                "workers": workers,
                "elapsed_s": time.perf_counter() - start,
            },
        )
        print(
            f"TX4_AUDIT_CHECKPOINT campaign={campaign} chunk={chunk_index + 1}/{n_chunks} "
            f"rows={len(shard_rows)} elapsed_s={time.perf_counter() - start:.1f}",
            flush=True,
        )


def assemble(campaign: str) -> None:
    shard_dir = RESULTS / "TX4_AUDIT_EXACT_SHARDS" / campaign
    paths = sorted(shard_dir.glob("shard_*.csv"))
    if not paths:
        raise RuntimeError(f"No shards found for {campaign}")
    frames = [pd.read_csv(path) for path in paths]
    df = pd.concat(frames, ignore_index=True)
    expected_conditions = 4096 if campaign == "QMC" else 5000
    if len(df) != expected_conditions * 16 or df.condition_id.nunique() != expected_conditions:
        raise RuntimeError(f"{campaign}: incomplete assembly {len(df)} rows/{df.condition_id.nunique()} conditions")
    df.to_parquet(RESULTS / f"TX4_{campaign}_ALL16_EXACT.parquet", index=False)
    df.to_csv(RESULTS / f"TX4_{campaign}_ALL16_EXACT.csv.gz", index=False, compression="gzip")
    print(json.dumps({"campaign": campaign, "rows": len(df), "conditions": expected_conditions}, indent=2), flush=True)


def main() -> None:
    requested = int(os.environ.get("TX4_AUDIT_WORKERS", "12"))
    workers = max(1, min(requested, (os.cpu_count() or 2) - 2))
    print(f"TX4_AUDIT_EXACT_START workers={workers} chunk_size={CHUNK_SIZE}", flush=True)
    for campaign in ("QMC", "MC"):
        run_campaign(campaign, workers)
        assemble(campaign)
    atomic_json(
        RESULTS / "TX4_AUDIT_EXACT_PROGRESS.json",
        {
            "status": "COMPLETE",
            "campaigns": {"QMC": {"conditions": 4096, "rows": 65536}, "MC": {"conditions": 5000, "rows": 80000}},
            "total_rows": 145536,
            "chunk_size": CHUNK_SIZE,
            "coordinates": "read from retained TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet; no samples regenerated",
        },
    )


if __name__ == "__main__":
    main()
