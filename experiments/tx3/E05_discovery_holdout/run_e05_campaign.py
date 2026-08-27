from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05 = ROOT / "experiments" / "tx3" / "E05_discovery_holdout"
sys.path[:0] = [str(ROOT), str(E03), str(E05)]

from campaign_core import action_vector, evaluate_operator, operating_points  # noqa: E402
from mechanism_consequence import (  # noqa: E402
    ModeSet,
    analyze_operating_point,
    modal_analysis,
    reference_catalog,
    state_names,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05_mechanism_consequence"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def worker(task: dict[str, Any]) -> dict[str, Any]:
    return analyze_operating_point(
        task["operating_point"],
        task["candidates"],
        task["reference_modes"],
        task["reference_records"],
        reference_band_hz=task["reference_band_hz"],
        tracking_band_hz=task["tracking_band_hz"],
        density_frequencies_hz=task["density_frequencies_hz"],
        sigma=task["sigma"],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", choices=("discovery", "holdout"), required=True)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    candidate_freeze = load_json(ARTIFACT / "preregistration" / "E05_CANDIDATE_FREEZE.json")
    numerical = load_json(ARTIFACT / "preregistration" / "E05_NUMERICAL_FREEZE.json")
    if args.split == "holdout" and not (ARTIFACT / "preregistration" / "E05_MODE_FAMILY_FREEZE.json").is_file():
        raise RuntimeError("holdout is locked until E05_MODE_FAMILY_FREEZE.json exists")
    candidates = {
        record["candidate_id"]: tuple(int(value[1:]) for value in record["actions"])
        for record in candidate_freeze["candidates"]
    }
    op00 = next(record for record in operating_points() if record["operating_point_id"] == "OP00")
    reference_point = evaluate_operator(op00, action_vector())
    if reference_point.status != "SUCCESS":
        raise RuntimeError(f"OP00 reference failed: {reference_point.status}")
    reference_modes: ModeSet = modal_analysis(reference_point)
    reference_records = reference_catalog(
        reference_modes,
        state_names(),
        band_hz=tuple(float(value) for value in numerical["reference_mode_band_hz"]),
    )
    catalog_path = ARTIFACT / "discovery" / "REFERENCE_MODE_CATALOG.json"
    catalog_payload = {
        "reference_operating_point": "OP00",
        "mode_count": len(reference_records),
        "modes": reference_records,
    }
    if args.split == "discovery":
        catalog_path.write_text(json.dumps(catalog_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        frozen_catalog = load_json(catalog_path)
        if len(frozen_catalog["modes"]) != len(reference_records):
            raise RuntimeError("OP00 reference mode catalog changed before holdout")
        for frozen, current in zip(frozen_catalog["modes"], reference_records, strict=True):
            if frozen["reference_mode_id"] != current["reference_mode_id"] or abs(frozen["frequency_hz"] - current["frequency_hz"]) > 1e-9:
                raise RuntimeError("OP00 reference mode identity changed before holdout")

    density_frequencies = np.geomspace(
        float(numerical["density_frequency_band_hz"][0]),
        float(numerical["density_frequency_band_hz"][1]),
        int(numerical["density_frequency_nodes"]),
    )
    selected_ops = [record for record in operating_points() if record["split"] == args.split]
    tasks = [
        {
            "operating_point": op,
            "candidates": candidates,
            "reference_modes": reference_modes,
            "reference_records": reference_records,
            "reference_band_hz": tuple(float(value) for value in numerical["reference_mode_band_hz"]),
            "tracking_band_hz": tuple(float(value) for value in numerical["tracking_target_band_hz"]),
            "density_frequencies_hz": density_frequencies,
            "sigma": float(numerical["sigma"]),
        }
        for op in selected_ops
    ]
    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(worker, task): task["operating_point"]["operating_point_id"] for task in tasks}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(f"{args.split}: {count}/{len(tasks)} {result['operating_point_id']} {result['status']}", flush=True)
    failures = [result for result in results if result["status"] != "SUCCESS"]
    if failures:
        (ARTIFACT / "logs" / f"{args.split}_failures.json").write_text(
            json.dumps(failures, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        raise RuntimeError(f"{len(failures)} E05 operating points failed")

    def frame(key: str) -> pd.DataFrame:
        rows = [row for result in results for row in result[key]]
        return pd.DataFrame(rows)

    tracking = frame("tracking_rows").sort_values(
        ["operating_point_id", "candidate_id", "reference_mode_id", "subset_order", "subset_key"]
    )
    consequences = frame("consequence_rows").sort_values(
        ["operating_point_id", "candidate_id", "reference_mode_id"]
    )
    system = frame("system_rows").sort_values(["operating_point_id", "subset_key"])
    density = frame("density_rows").sort_values(["operating_point_id", "candidate_id", "grid_index"])
    output = ARTIFACT / args.split
    tracking.to_parquet(output / "mode_tracking_vertices.parquet", index=False)
    consequences.to_parquet(output / "mode_consequences_all_families.parquet", index=False)
    system.to_parquet(output / "system_vertex_metrics.parquet", index=False)
    density.to_parquet(ARTIFACT / "density" / f"{args.split}_externality_density.parquet", index=False)
    compute = {
        "split": args.split,
        "operating_points": len(selected_ops),
        "workers": args.workers,
        "wall_time_s": time.perf_counter() - started,
        "operator_builds": sum(int(result["operator_builds"]) for result in results),
        "reference_modes": len(reference_records),
        "tracking_rows": len(tracking),
        "consequence_rows": len(consequences),
        "density_rows": len(density),
        "status": "SUCCESS",
    }
    (ARTIFACT / "logs" / f"{args.split}_compute.json").write_text(
        json.dumps(compute, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(compute, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
