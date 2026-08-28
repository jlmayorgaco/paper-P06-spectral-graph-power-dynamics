from __future__ import annotations

import argparse
import json
import math
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
import psutil


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E03), str(E05D)]

from campaign_core import action_vector, coalition_externality  # noqa: E402
from weak_grid import (  # noqa: E402
    GeneralizedModes,
    complex_externality,
    damping_ratio,
    evaluate_weak_grid_operator,
    generalized_modes,
    maximum_physical_real_part,
    mobius_vertices,
    parse_subset_key,
    stable_lower_order_prediction,
    subset_key,
    track_one_hungarian,
    union_vertices,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
PREREG = ARTIFACT / "preregistration"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def worker(task: dict[str, Any]) -> dict[str, Any]:
    seed = task["seed"]
    kappa_end = float(task["kappa_end"])
    tau_grid = [float(value) for value in task["stress"]["tau_grid"]]
    candidates = task["candidates"]
    frozen_empty = task["frozen_empty"]
    bmac_min = float(task["numerical"]["stress_and_action_tracking_bmac_minimum"])
    residual_max = float(task["numerical"]["eigenpair_normalized_residual_maximum"])
    screen = float(task["gates"]["screen_damping_ratio"])
    started = time.perf_counter()
    process = psutil.Process()
    peak_rss = process.memory_info().rss
    cache: dict[tuple[float, tuple[int, ...], float], tuple[Any, GeneralizedModes | None, str]] = {}
    build_rows: list[dict[str, Any]] = []
    cache_requests = 0

    def cache_key(kappa: float, subset: tuple[int, ...], beta: float) -> tuple[float, tuple[int, ...], float]:
        if beta == 0.0 or not subset:
            return (round(float(kappa), 12), (), 0.0)
        return (round(float(kappa), 12), tuple(subset), round(float(beta), 12))

    def get(kappa: float, subset: tuple[int, ...] = (), beta: float = 0.0) -> tuple[Any, GeneralizedModes | None, str]:
        nonlocal cache_requests, peak_rss
        cache_requests += 1
        key = cache_key(kappa, subset, beta)
        if key in cache:
            return cache[key]
        alpha = action_vector({f"A{value}": float(beta) for value in subset})
        point = evaluate_weak_grid_operator(seed, kappa, alpha)
        modes = None
        status = point.status
        eig_exception = None
        if point.status == "SUCCESS":
            try:
                modes = generalized_modes(point)
                if float(np.max(modes.residuals)) > residual_max:
                    status = "EIGENPAIR_RESIDUAL_FAIL"
            except Exception as exc:
                status = "EIG_FAIL"
                eig_exception = repr(exc)
        cache[key] = (point, modes, status)
        build_rows.append({
            "seed_id": seed["seed_id"],
            "tau": None,
            "kappa_grid": float(kappa),
            "coalition": "UNION",
            "order": len(subset),
            "subset": subset_key(subset),
            "beta": float(beta),
            "status": status,
            "eig_exception": eig_exception,
            "maximum_eigenpair_residual": float(np.max(modes.residuals)) if modes is not None else None,
            "maximum_physical_real_part_per_s": maximum_physical_real_part(modes) if modes is not None else None,
            **point.metadata,
        })
        peak_rss = max(peak_rss, process.memory_info().rss)
        return cache[key]

    def track_action(
        kappa: float, empty_modes: GeneralizedModes, empty_index: int, subset: tuple[int, ...]
    ) -> tuple[int, float, int, str]:
        if not subset:
            return empty_index, 1.0, 1, "SUCCESS"
        last_status = "MODE_TRACK_FAIL"
        for count in (5, 9, 17):
            previous_modes = empty_modes
            previous_index = empty_index
            minimum = 1.0
            success = True
            for beta in np.linspace(0.0, 1.0, count)[1:]:
                _, modes, status = get(kappa, subset, float(beta))
                if status != "SUCCESS" or modes is None:
                    success = False
                    last_status = "NONSMOOTH_REGIME" if status == "LIMITER_ACTIVE" else status
                    break
                match = track_one_hungarian(previous_modes, modes, previous_index)
                if match is None or match[1] < bmac_min:
                    success = False
                    last_status = "MODE_TRACK_FAIL"
                    break
                previous_index, score = match
                previous_modes = modes
                minimum = min(minimum, score)
            if success:
                return previous_index, minimum, count, "SUCCESS"
        return -1, 0.0, 17, last_status

    coalitions = [parse_subset_key(record["coalition_key"]) for record in candidates]
    subsets = union_vertices(coalitions)
    subset_rows: list[dict[str, Any]] = []
    externality_rows: list[dict[str, Any]] = []
    for tau in tau_grid:
        kappa = 1.0 + tau * (kappa_end - 1.0)
        point0, modes0, status0 = get(kappa)
        freeze = frozen_empty[round(tau, 6)]
        if status0 != "SUCCESS" or modes0 is None:
            for candidate in candidates:
                externality_rows.append({
                    "seed_id": seed["seed_id"], "tau": tau, "kappa_grid": kappa,
                    "coalition": candidate["coalition_key"], "coalition_key": candidate["coalition_key"],
                    "candidate_id": candidate["candidate_id"], "order": candidate["order"], "status": status0,
                })
            continue
        empty_index = int(freeze["mode_index"])
        if empty_index >= len(modes0.eigenvalues):
            raise RuntimeError(f"{seed['seed_id']} frozen empty index out of range at tau={tau}")
        empty_value = complex(modes0.eigenvalues[empty_index])
        frozen_value = complex(float(freeze["lambda_real"]), float(freeze["lambda_imag"]))
        if abs(empty_value - frozen_value) > 1e-8:
            raise RuntimeError(f"{seed['seed_id']} frozen empty eigenvalue mismatch at tau={tau}")
        tracked_subsets: dict[tuple[int, ...], dict[str, Any]] = {}
        for subset in subsets:
            target_index, action_bmac, homotopy_nodes, status = track_action(kappa, modes0, empty_index, subset)
            if status != "SUCCESS" or target_index < 0:
                tracked_subsets[subset] = {"status": status}
                subset_rows.append({
                    "seed_id": seed["seed_id"], "tau": tau, "kappa_grid": kappa,
                    "min_SCR": float(point0.metadata["min_SCR"]), "SCR36": float(point0.metadata["SCR36"]),
                    "SCR37": float(point0.metadata["SCR37"]), "SCR38": float(point0.metadata["SCR38"]),
                    "coalition": "UNION", "order": len(subset), "subset": subset_key(subset),
                    "BMAC": 0.0, "status": status,
                })
                continue
            target_point, target_modes, endpoint_status = get(kappa, subset, 1.0 if subset else 0.0)
            if endpoint_status != "SUCCESS" or target_modes is None:
                tracked_subsets[subset] = {"status": endpoint_status}
                continue
            value = complex(target_modes.eigenvalues[target_index])
            combined_bmac = min(float(freeze["BMAC"]), float(action_bmac))
            record = {
                "seed_id": seed["seed_id"], "tau": tau, "kappa_grid": kappa,
                "min_SCR": float(point0.metadata["min_SCR"]), "SCR36": float(point0.metadata["SCR36"]),
                "SCR37": float(point0.metadata["SCR37"]), "SCR38": float(point0.metadata["SCR38"]),
                "coalition": "UNION", "order": len(subset), "subset": subset_key(subset),
                "lambda_real": float(value.real), "lambda_imag": float(value.imag),
                "frequency_hz": float(value.imag / (2.0 * math.pi)), "damping_ratio": damping_ratio(value),
                "BMAC": combined_bmac, "stress_BMAC": float(freeze["BMAC"]), "action_BMAC": float(action_bmac),
                "stress_refinement_depth": int(freeze["stress_refinement_depth"]),
                "action_homotopy_nodes": homotopy_nodes,
                "eigenpair_residual": float(target_modes.residuals[target_index]),
                "kappa_lambda": float(target_modes.condition_numbers[target_index]),
                "nearest_modal_separation_per_s": float(target_modes.nearest_separations[target_index]),
                "status": "SUCCESS",
            }
            tracked_subsets[subset] = {**record, "modes": target_modes, "mode_index": target_index, "point": target_point}
            subset_rows.append(record)
        for candidate in candidates:
            coalition = parse_subset_key(candidate["coalition_key"])
            required = mobius_vertices(coalition)
            base = {
                "seed_id": seed["seed_id"], "tau": tau, "kappa_grid": kappa,
                "min_SCR": float(point0.metadata["min_SCR"]), "SCR36": float(point0.metadata["SCR36"]),
                "SCR37": float(point0.metadata["SCR37"]), "SCR38": float(point0.metadata["SCR38"]),
                "candidate_id": candidate["candidate_id"], "coalition": candidate["coalition_key"],
                "coalition_key": candidate["coalition_key"], "order": candidate["order"], "subset": "MOBIUS",
            }
            failed = next((tracked_subsets[subset]["status"] for subset in required if tracked_subsets[subset]["status"] != "SUCCESS"), None)
            if failed is not None:
                externality_rows.append({**base, "status": failed})
                continue
            values = {
                subset: complex(tracked_subsets[subset]["lambda_real"], tracked_subsets[subset]["lambda_imag"])
                for subset in required
            }
            damping = {subset: float(tracked_subsets[subset]["damping_ratio"]) for subset in required}
            delta = complex_externality(values, coalition)
            delta_zeta = coalition_externality(damping, coalition)
            lower = stable_lower_order_prediction(values, coalition)
            full = values[tuple(sorted(coalition))]
            lambda0 = values[()]
            closure = full - (lower + delta)
            rho_comp = delta.real / (-lower.real) if lower.real < 0.0 else math.nan
            rho_baseline = delta.real / (-lambda0.real) if lambda0.real < 0.0 else math.nan
            zeta_lower = damping_ratio(lower)
            zeta_full = damping_ratio(full)
            hard_failure = bool(lower.real < 0.0 and full.real >= 0.0)
            screen_failure = bool(zeta_lower >= screen and zeta_full < screen)
            full_record = tracked_subsets[tuple(sorted(coalition))]
            status = "LOWER_ORDER_UNSTABLE" if lower.real >= 0.0 else "SUCCESS"
            externality_rows.append({
                **base,
                "lambda_real": float(full.real), "lambda_imag": float(full.imag),
                "frequency_hz": float(full.imag / (2.0 * math.pi)), "damping_ratio": zeta_full,
                "BMAC": float(min(tracked_subsets[subset]["BMAC"] for subset in required)),
                "eigenpair_residual": float(full_record["eigenpair_residual"]),
                "kappa_lambda": float(full_record["kappa_lambda"]),
                "nearest_modal_separation_per_s": float(full_record["nearest_modal_separation_per_s"]),
                "critical_lambda_real": float(lambda0.real), "critical_lambda_imag": float(lambda0.imag),
                "critical_damping_ratio": damping_ratio(lambda0),
                "critical_kappa_lambda": float(tracked_subsets[()]["kappa_lambda"]),
                "critical_nearest_modal_separation_per_s": float(tracked_subsets[()]["nearest_modal_separation_per_s"]),
                "Delta_lambda_real": float(delta.real), "Delta_lambda_imag": float(delta.imag),
                "Delta_lambda_absolute": float(abs(delta)), "Delta_zeta": float(delta_zeta),
                "lambda_lower_real": float(lower.real), "lambda_lower_imag": float(lower.imag),
                "zeta_lower": float(zeta_lower), "rho_comp": float(rho_comp),
                "rho_baseline": float(rho_baseline), "hard_failure": hard_failure,
                "screen_failure": screen_failure,
                "composition_closure_real": float(closure.real),
                "composition_closure_imag": float(closure.imag),
                "destabilizing_material_point": bool(rho_comp >= 0.25 and delta_zeta <= -0.0025),
                "status": status,
            })
    return {
        "seed_id": seed["seed_id"], "subsets": subset_rows, "externalities": externality_rows,
        "builds": build_rows, "wall_time_s": time.perf_counter() - started,
        "physical_build_count": len(cache), "eigensolve_count": sum(row["maximum_eigenpair_residual"] is not None for row in build_rows),
        "cache_requests": cache_requests, "cache_hits": cache_requests - len(cache), "peak_rss_bytes": peak_rss,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    endpoint_freeze = PREREG / "E05D_STRESS_ENDPOINT_FREEZE.json"
    mode_freeze = PREREG / "E05D_CRITICAL_MODE_FREEZE.json"
    if not endpoint_freeze.is_file() or not mode_freeze.is_file():
        raise FileNotFoundError("frozen baseline endpoints and critical modes are required before coalition outcomes")
    candidates = load_json(PREREG / "E05D_CANDIDATE_FREEZE.json")["candidates"]
    seeds = load_json(PREREG / "E05D_HOLDOUT_FREEZE.json")["operating_points"]
    stress = load_json(PREREG / "E05D_STRESS_RULE_FREEZE.json")
    numerical = load_json(PREREG / "E05D_NUMERICAL_FREEZE.json")
    gates = load_json(PREREG / "E05D_C3C_GATE_FREEZE.json")
    limits = pd.read_parquet(ARTIFACT / "baseline_calibration" / "E05D_BASELINE_STRESS_LIMITS.parquet").set_index("seed_id")
    frozen = pd.read_parquet(PREREG / "e05d_empty_mode_tracking_freeze.parquet")
    tasks = []
    for seed in seeds:
        rows = frozen[frozen["seed_id"] == seed["seed_id"]]
        tasks.append({
            "seed": seed, "kappa_end": float(limits.loc[seed["seed_id"], "kappa_end"]),
            "candidates": candidates, "stress": stress, "numerical": numerical, "gates": gates,
            "frozen_empty": {round(float(row.tau), 6): row._asdict() for row in rows.itertuples(index=False)},
        })
    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(worker, task): task["seed"]["seed_id"] for task in tasks}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(f"{count}/8 {result['seed_id']} {result['physical_build_count']} builds {result['wall_time_s']:.1f}s", flush=True)
    results.sort(key=lambda value: value["seed_id"])
    subsets = pd.DataFrame([row for result in results for row in result["subsets"]])
    externalities = pd.DataFrame([row for result in results for row in result["externalities"]])
    builds = pd.DataFrame([row for result in results for row in result["builds"]])
    tau_lookup = {(round(float(row.kappa_grid), 12), row.seed_id): float(row.tau) for row in subsets.itertuples(index=False)}
    builds["tau"] = [tau_lookup.get((round(float(kappa), 12), seed)) for kappa, seed in zip(builds["kappa_grid"], builds["seed_id"], strict=True)]
    output = ARTIFACT / "holdout"
    output.mkdir(parents=True, exist_ok=True)
    subsets.to_parquet(output / "e05d_subset_poles.parquet", index=False)
    subsets.to_parquet(output / "e05d_mode_tracking.parquet", index=False)
    externalities.to_parquet(output / "e05d_connected_externalities.parquet", index=False)
    externalities.to_parquet(output / "e05d_margin_consumption.parquet", index=False)
    externalities[[
        "seed_id", "tau", "kappa_grid", "min_SCR", "SCR36", "SCR37", "SCR38",
        "candidate_id", "coalition", "order", "lambda_lower_real", "lambda_lower_imag",
        "lambda_real", "lambda_imag", "zeta_lower", "damping_ratio", "hard_failure", "screen_failure", "status",
    ]].to_parquet(output / "e05d_composition_failures.parquet", index=False)
    externalities[[
        "seed_id", "tau", "kappa_grid", "min_SCR", "candidate_id", "coalition", "order",
        "rho_comp", "rho_baseline", "critical_kappa_lambda", "kappa_lambda",
        "critical_nearest_modal_separation_per_s", "nearest_modal_separation_per_s", "status",
    ]].to_parquet(output / "e05d_conditioning.parquet", index=False)
    builds.to_parquet(output / "e05d_failure_ledger.parquet", index=False)
    total_requests = sum(result["cache_requests"] for result in results)
    total_hits = sum(result["cache_hits"] for result in results)
    compute = {
        "seed_count": 8, "tau_count": 9, "candidate_count": 15,
        "unique_endpoint_subset_count": len(union_vertices([parse_subset_key(record["coalition_key"]) for record in candidates])),
        "physical_build_count": sum(result["physical_build_count"] for result in results),
        "pf_count": sum(result["physical_build_count"] for result in results),
        "eigensolve_count": sum(result["eigensolve_count"] for result in results),
        "cache_requests": total_requests, "cache_hits": total_hits,
        "cache_hit_rate": total_hits / total_requests,
        "externality_rows": len(externalities), "successful_externality_rows": int(externalities["status"].eq("SUCCESS").sum()),
        "subset_pole_rows": len(subsets), "build_status_counts": builds["status"].value_counts(dropna=False).to_dict(),
        "peak_worker_rss_bytes": max(result["peak_rss_bytes"] for result in results),
        "workers": args.workers, "wall_time_s": time.perf_counter() - started,
    }
    (ARTIFACT / "logs").mkdir(parents=True, exist_ok=True)
    (ARTIFACT / "logs" / "e05d_campaign_compute.json").write_text(
        json.dumps(compute, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(compute, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
