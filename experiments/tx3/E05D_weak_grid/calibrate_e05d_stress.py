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
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E05D)]

from weak_grid import (  # noqa: E402
    baseline_boundary,
    critical_mode_index,
    damping_ratio,
    evaluate_weak_grid_operator,
    generalized_modes,
    maximum_physical_real_part,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
PREREG = ARTIFACT / "preregistration"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def evaluate(seed: dict[str, Any], kappa: float) -> dict[str, Any]:
    point = evaluate_weak_grid_operator(seed, kappa, {})
    modes = None
    eig_exception = None
    physical_status = point.status
    if point.status == "SUCCESS":
        try:
            modes = generalized_modes(point)
        except Exception as exc:
            physical_status = "EIG_FAIL"
            eig_exception = repr(exc)
    reason = baseline_boundary(point, modes)
    critical_index = None
    critical_value = None
    if modes is not None:
        try:
            critical_index = critical_mode_index(modes)
            critical_value = complex(modes.eigenvalues[critical_index])
        except ValueError:
            pass
    return {
        "seed_id": seed["seed_id"],
        "kappa_grid": float(kappa),
        "status": physical_status,
        "boundary_reason": reason,
        "boundary_reached": reason is not None,
        "critical_mode_index_at_point": critical_index,
        "minimum_oscillatory_damping_ratio": damping_ratio(critical_value) if critical_value is not None else None,
        "critical_lambda_real_per_s": float(critical_value.real) if critical_value is not None else None,
        "critical_lambda_imag_per_s": float(critical_value.imag) if critical_value is not None else None,
        "maximum_physical_real_part_per_s": maximum_physical_real_part(modes) if modes is not None else None,
        "maximum_eigenpair_residual": float(np.max(modes.residuals)) if modes is not None else None,
        "eig_exception": eig_exception,
        **point.metadata,
    }


def worker(task: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    seed, stress = task
    coarse = [float(value) for value in stress["coarse_kappa_grid"]]
    tolerance = float(stress["boundary_bisection_tolerance"])
    trace: list[dict[str, Any]] = []
    first = evaluate(seed, 1.0)
    trace.append({**first, "calibration_phase": "coarse"})
    if first["boundary_reached"]:
        return {"seed_id": seed["seed_id"], "status": "INVALID_AT_KAPPA_1", "trace": trace}
    low = 1.0
    first_failure = None
    for kappa in coarse[1:]:
        row = evaluate(seed, kappa)
        trace.append({**row, "calibration_phase": "coarse"})
        if row["boundary_reached"]:
            first_failure = row
            break
        low = kappa
    if first_failure is None:
        kappa_limit = float(stress["kappa_cap"])
        kappa_end = kappa_limit
        endpoint_reason = "CAP_REACHED"
        endpoint_classification = "CAP_REACHED_WITHOUT_REGISTERED_BOUNDARY"
        safe_fraction = 1.0
    else:
        high = float(first_failure["kappa_grid"])
        while high - low > tolerance:
            midpoint = 0.5 * (low + high)
            row = evaluate(seed, midpoint)
            trace.append({**row, "calibration_phase": "bisection"})
            if row["boundary_reached"]:
                high = midpoint
                first_failure = row
            else:
                low = midpoint
        kappa_limit = high
        endpoint_reason = str(first_failure["boundary_reason"])
        if endpoint_reason in {"DAMPING_SCREEN", "BASELINE_UNSTABLE"}:
            safe_fraction = 0.90
            endpoint_classification = "DYNAMIC_MARGIN_BOUNDARY"
        else:
            safe_fraction = 0.95
            endpoint_classification = "PHYSICAL_LIMIT_PRECEDES_DYNAMIC_MARGIN"
        kappa_end = 1.0 + safe_fraction * (kappa_limit - 1.0)
    endpoint = evaluate(seed, kappa_end)
    trace.append({**endpoint, "calibration_phase": "safe_endpoint_check"})
    return {
        "seed_id": seed["seed_id"],
        "status": "SUCCESS" if not endpoint["boundary_reached"] else "NUMERICAL_UNRESOLVED",
        "kappa_limit": float(kappa_limit),
        "kappa_end": float(kappa_end),
        "endpoint_reason": endpoint_reason,
        "endpoint_classification": endpoint_classification,
        "safe_endpoint_fraction": float(safe_fraction),
        "baseline_critical_damping_ratio": float(first["minimum_oscillatory_damping_ratio"]),
        "endpoint_critical_damping_ratio": float(endpoint["minimum_oscillatory_damping_ratio"]),
        "baseline_min_SCR": float(first["min_SCR"]),
        "endpoint_min_SCR": float(endpoint["min_SCR"]),
        "endpoint_SCR36": float(endpoint["SCR36"]),
        "endpoint_SCR37": float(endpoint["SCR37"]),
        "endpoint_SCR38": float(endpoint["SCR38"]),
        "trace": trace,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    if (ARTIFACT / "holdout" / "e05d_connected_externalities.parquet").exists():
        raise RuntimeError("coalition outcomes already exist; baseline calibration cannot be rerun")
    seeds = load_json(PREREG / "E05D_HOLDOUT_FREEZE.json")["operating_points"]
    stress = load_json(PREREG / "E05D_STRESS_RULE_FREEZE.json")
    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(worker, (seed, stress)) for seed in seeds]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(result["seed_id"], result["status"], result.get("kappa_end"), flush=True)
    results.sort(key=lambda value: value["seed_id"])
    if any(result["status"] != "SUCCESS" for result in results):
        raise RuntimeError(f"baseline-only calibration failed: {[(r['seed_id'], r['status']) for r in results]}")
    limits = pd.DataFrame([{key: value for key, value in result.items() if key != "trace"} for result in results])
    trace = pd.DataFrame([row for result in results for row in result["trace"]])
    output = ARTIFACT / "baseline_calibration"
    output.mkdir(parents=True, exist_ok=True)
    (ARTIFACT / "tables").mkdir(parents=True, exist_ok=True)
    limits.to_parquet(output / "E05D_BASELINE_STRESS_LIMITS.parquet", index=False)
    limits.to_parquet(ARTIFACT / "tables" / "e05d_seed_stress_limits.parquet", index=False)
    trace.to_parquet(output / "E05D_BASELINE_CALIBRATION_TRACE.parquet", index=False)
    payload = {
        "status": "SUCCESS",
        "empty_coalition_only": True,
        "coalition_evaluations": 0,
        "seed_count": len(results),
        "kappa_end_range": [float(limits["kappa_end"].min()), float(limits["kappa_end"].max())],
        "endpoint_reasons": limits["endpoint_reason"].value_counts().to_dict(),
        "evaluation_count": len(trace),
        "wall_time_s": time.perf_counter() - started,
    }
    (ARTIFACT / "logs").mkdir(parents=True, exist_ok=True)
    (ARTIFACT / "logs" / "baseline_calibration_compute.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
