from __future__ import annotations

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

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
sys.path[:0] = [str(ROOT), str(E05C)]

from margin_conditioned import (  # noqa: E402
    evaluate_stressed_operator,
    generalized_modes,
    hard_stop_reason,
    maximum_physical_real_part,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05C_margin_conditioned"
PREREG = ARTIFACT / "preregistration"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def evaluate(seed: dict[str, Any], gamma: float) -> dict[str, Any]:
    point = evaluate_stressed_operator(seed, gamma, {})
    modes = None
    eig_exception = None
    if point.status == "SUCCESS":
        try:
            modes = generalized_modes(point)
        except Exception as exc:
            eig_exception = repr(exc)
    reason = hard_stop_reason(point, modes)
    return {
        "seed_id": seed["seed_id"],
        "gamma": float(gamma),
        "status": "SUCCESS" if reason is None else reason,
        "hard_stop": reason is not None,
        "maximum_physical_real_part_per_s": maximum_physical_real_part(modes) if modes is not None else None,
        "maximum_eigenpair_residual": float(modes.residuals.max()) if modes is not None else None,
        "eig_exception": eig_exception,
        **point.metadata,
    }


def worker(task: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    seed, stress = task
    increment = float(stress["coarse_increment"])
    cap = float(stress["gamma_cap"])
    tolerance = float(stress["boundary_bisection_tolerance"])
    trace = []
    gamma = 1.0
    previous_success = None
    first_failure = None
    while gamma <= cap + 1e-12:
        row = evaluate(seed, round(gamma, 12))
        trace.append(row)
        if row["hard_stop"]:
            first_failure = row
            break
        previous_success = row
        gamma += increment
    if previous_success is None:
        return {"seed_id": seed["seed_id"], "trace": trace, "status": "INVALID_AT_GAMMA_1"}
    if first_failure is None:
        gamma_limit = cap
        gamma_end = cap
        endpoint_reason = "NO_HARD_STOP_THROUGH_CAP"
        bracket_low = cap
        bracket_high = cap
    else:
        low = float(previous_success["gamma"])
        high = float(first_failure["gamma"])
        while high - low > tolerance:
            midpoint = 0.5 * (low + high)
            row = evaluate(seed, midpoint)
            trace.append(row)
            if row["hard_stop"]:
                high = midpoint
            else:
                low = midpoint
        gamma_limit = high
        gamma_end = 1.0 + float(stress["safe_endpoint_fraction"]) * (gamma_limit - 1.0)
        endpoint_reason = str(first_failure["status"])
        bracket_low, bracket_high = low, high
    endpoint = evaluate(seed, gamma_end)
    trace.append({**endpoint, "endpoint_check": True})
    status = "SUCCESS" if not endpoint["hard_stop"] else "SAFE_ENDPOINT_FAIL"
    return {
        "seed_id": seed["seed_id"],
        "status": status,
        "gamma_limit": float(gamma_limit),
        "gamma_end": float(gamma_end),
        "endpoint_reason": endpoint_reason,
        "boundary_bracket_low": float(bracket_low),
        "boundary_bracket_high": float(bracket_high),
        "endpoint_voltage_min_pu": endpoint.get("voltage_min_pu"),
        "endpoint_voltage_max_pu": endpoint.get("voltage_max_pu"),
        "endpoint_maximum_physical_real_part_per_s": endpoint.get("maximum_physical_real_part_per_s"),
        "trace": trace,
    }


def main() -> int:
    if (ARTIFACT / "holdout" / "e05c_connected_pole_externalities.parquet").exists():
        raise RuntimeError("coalition results already exist; baseline calibration cannot be rerun")
    seeds = load_json(PREREG / "E05C_BASE_HOLDOUT_FREEZE.json")["operating_points"]
    stress = load_json(PREREG / "E05C_STRESS_RULE_FREEZE.json")
    started = time.perf_counter()
    results = []
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(worker, (seed, stress)): seed["seed_id"] for seed in seeds}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(result["seed_id"], result["status"], result.get("gamma_end"), flush=True)
    results.sort(key=lambda value: value["seed_id"])
    if any(result["status"] != "SUCCESS" for result in results):
        raise RuntimeError(f"baseline-only calibration failed: {[(r['seed_id'], r['status']) for r in results]}")
    limits = pd.DataFrame([{key: value for key, value in result.items() if key != "trace"} for result in results])
    trace = pd.DataFrame([row for result in results for row in result["trace"]])
    limits.to_parquet(ARTIFACT / "baseline_calibration" / "E05C_BASELINE_STRESS_LIMITS.parquet", index=False)
    limits.to_parquet(ARTIFACT / "tables" / "e05c_seed_stress_limits.parquet", index=False)
    trace.to_parquet(ARTIFACT / "baseline_calibration" / "E05C_BASELINE_CALIBRATION_TRACE.parquet", index=False)
    payload = {
        "status": "SUCCESS",
        "baseline_only": True,
        "coalition_vertices_evaluated": 0,
        "seed_count": len(results),
        "gamma_end_minimum": float(limits["gamma_end"].min()),
        "gamma_end_maximum": float(limits["gamma_end"].max()),
        "endpoint_reasons": limits["endpoint_reason"].value_counts().to_dict(),
        "wall_time_s": time.perf_counter() - started,
    }
    (ARTIFACT / "logs" / "baseline_calibration_compute.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
