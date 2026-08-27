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


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
sys.path[:0] = [str(ROOT), str(E03), str(E05C)]

from campaign_core import action_vector, coalition_externality  # noqa: E402
from margin_conditioned import (  # noqa: E402
    GeneralizedModes,
    complex_externality,
    damping_ratio,
    evaluate_stressed_operator,
    generalized_modes,
    maximum_physical_real_part,
    mobius_vertices,
    parse_subset_key,
    pole_factor,
    stable_lower_order_prediction,
    subset_key,
    track_one_hungarian,
    union_vertices,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05C_margin_conditioned"
PREREG = ARTIFACT / "preregistration"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def reference_modes() -> GeneralizedModes:
    seed = {
        "seed_id": "OP00",
        "load_scale": 1.0,
        "gfl_dispatch_scale": 1.0,
        "reactive_load_scale": 1.0,
        "redispatch_mw": 0.0,
    }
    point = evaluate_stressed_operator(seed, 1.0, {})
    if point.status != "SUCCESS":
        raise RuntimeError(f"OP00 reference failed: {point.status}")
    return generalized_modes(point)


def worker(task: dict[str, Any]) -> dict[str, Any]:
    seed = task["seed"]
    gamma_end = float(task["gamma_end"])
    tau_grid = [float(value) for value in task["stress"]["tau_grid"]]
    candidates = task["candidates"]
    locks = task["locks"]
    development_modes: GeneralizedModes = task["development_modes"]
    bmac_min = float(task["numerical"]["stress_and_action_tracking_bmac_minimum"])
    residual_max = float(task["numerical"]["eigenpair_normalized_residual_maximum"])
    screen = float(task["gates"]["screen_damping_ratio"])
    started = time.perf_counter()
    cache: dict[tuple[float, tuple[int, ...], float], tuple[Any, GeneralizedModes | None]] = {}
    build_rows: list[dict[str, Any]] = []
    pole_cache: dict[tuple[float, tuple[int, ...], float], float] = {}

    def cache_key(gamma: float, subset: tuple[int, ...], beta: float) -> tuple[float, tuple[int, ...], float]:
        if beta == 0.0 or not subset:
            return (round(float(gamma), 12), (), 0.0)
        return (round(float(gamma), 12), tuple(subset), round(float(beta), 12))

    def get(gamma: float, subset: tuple[int, ...] = (), beta: float = 0.0) -> tuple[Any, GeneralizedModes | None]:
        key = cache_key(gamma, subset, beta)
        if key in cache:
            return cache[key]
        alpha = action_vector({f"A{value}": float(beta) for value in subset})
        point = evaluate_stressed_operator(seed, gamma, alpha)
        modes = None
        eig_exception = None
        status = point.status
        if point.status == "SUCCESS":
            try:
                modes = generalized_modes(point)
                if float(np.max(modes.residuals)) > residual_max:
                    status = "EIGENPAIR_RESIDUAL_FAIL"
            except Exception as exc:
                status = "EIG_FAIL"
                eig_exception = repr(exc)
        cache[key] = (point, modes)
        build_rows.append(
            {
                "seed_id": seed["seed_id"],
                "gamma": float(gamma),
                "subset": subset_key(subset),
                "beta": float(beta),
                "status": status,
                "eig_exception": eig_exception,
                "maximum_eigenpair_residual": float(np.max(modes.residuals)) if modes is not None else None,
                "maximum_physical_real_part_per_s": maximum_physical_real_part(modes) if modes is not None else None,
                **point.metadata,
            }
        )
        return point, modes

    def get_pole(gamma: float, subset: tuple[int, ...]) -> float:
        key = cache_key(gamma, subset, 1.0 if subset else 0.0)
        if key not in pole_cache:
            _, modes = get(gamma, subset, 1.0 if subset else 0.0)
            if modes is None:
                raise ValueError("pole factor unavailable")
            pole_cache[key] = pole_factor(modes, int(task["numerical"]["dynamic_pole_factor_nodes"]))
        return pole_cache[key]

    def track_stress(
        source_modes: GeneralizedModes,
        source_index: int,
        source_gamma: float,
        target_gamma: float,
        depth: int = 0,
    ) -> tuple[int, float, int] | None:
        _, target_modes = get(target_gamma)
        if target_modes is None:
            return None
        match = track_one_hungarian(source_modes, target_modes, source_index)
        if match is not None and match[1] >= bmac_min:
            return match[0], match[1], depth
        if depth >= 2:
            return None
        midpoint = 0.5 * (source_gamma + target_gamma)
        _, middle_modes = get(midpoint)
        if middle_modes is None:
            return None
        first = track_one_hungarian(source_modes, middle_modes, source_index)
        if first is None or first[1] < bmac_min:
            return None
        second = track_stress(middle_modes, first[0], midpoint, target_gamma, depth + 1)
        if second is None:
            return None
        return second[0], min(first[1], second[1]), max(depth + 1, second[2])

    def track_action(
        gamma: float, empty_modes: GeneralizedModes, empty_index: int, subset: tuple[int, ...]
    ) -> tuple[int, float, int, str] | None:
        if not subset:
            return empty_index, 1.0, 1, "SUCCESS"
        for count in (5, 9, 17):
            nodes = np.linspace(0.0, 1.0, count)
            previous_modes = empty_modes
            previous_index = empty_index
            minimum = 1.0
            success = True
            failure_status = "MODE_TRACK_FAIL"
            for beta in nodes[1:]:
                point, modes = get(gamma, subset, float(beta))
                if point.status != "SUCCESS" or modes is None:
                    success = False
                    failure_status = point.status if point.status != "SUCCESS" else "EIG_FAIL"
                    break
                match = track_one_hungarian(previous_modes, modes, previous_index)
                if match is None or match[1] < bmac_min:
                    success = False
                    failure_status = "MODE_TRACK_FAIL"
                    break
                previous_index, score = match
                previous_modes = modes
                minimum = min(minimum, score)
            if success:
                return previous_index, minimum, count, "SUCCESS"
        return None if failure_status == "MODE_TRACK_FAIL" else (-1, 0.0, 17, failure_status)

    track_rows: list[dict[str, Any]] = []
    externality_rows: list[dict[str, Any]] = []
    previous_empty: dict[str, tuple[float, GeneralizedModes, int]] = {}
    for tau in tau_grid:
        gamma = 1.0 + tau * (gamma_end - 1.0)
        empty_point, empty_modes = get(gamma)
        if empty_point.status != "SUCCESS" or empty_modes is None:
            for candidate in candidates:
                externality_rows.append(
                    {
                        "seed_id": seed["seed_id"], "tau": tau, "gamma": gamma,
                        "candidate_id": candidate["candidate_id"], "coalition_key": candidate["coalition_key"],
                        "order": candidate["order"], "status": empty_point.status,
                    }
                )
            continue
        for candidate in candidates:
            lock = locks[candidate["candidate_id"]]
            base = {
                "seed_id": seed["seed_id"],
                "tau": tau,
                "gamma": gamma,
                "candidate_id": candidate["candidate_id"],
                "coalition_key": candidate["coalition_key"],
                "order": candidate["order"],
                "modal_family_id": lock.get("modal_family_id"),
            }
            if lock["lock_status"] != "LOCKED":
                externality_rows.append({**base, "status": "NO_LOCKABLE_MODE"})
                continue
            if tau == tau_grid[0]:
                match = track_one_hungarian(development_modes, empty_modes, int(lock["reference_mode_index"]))
                if match is None or match[1] < bmac_min:
                    externality_rows.append({**base, "status": "MODE_TRACK_FAIL"})
                    continue
                empty_index, stress_bmac, stress_depth = match[0], match[1], 0
            else:
                previous = previous_empty.get(candidate["candidate_id"])
                if previous is None:
                    externality_rows.append({**base, "status": "MODE_TRACK_FAIL"})
                    continue
                tracked = track_stress(previous[1], previous[2], previous[0], gamma)
                if tracked is None:
                    externality_rows.append({**base, "status": "MODE_TRACK_FAIL"})
                    continue
                empty_index, stress_bmac, stress_depth = tracked
            previous_empty[candidate["candidate_id"]] = (gamma, empty_modes, empty_index)
            coalition = parse_subset_key(candidate["coalition_key"])
            values: dict[tuple[int, ...], complex] = {}
            damping: dict[tuple[int, ...], float] = {}
            modes_for_vertex: dict[tuple[int, ...], GeneralizedModes] = {}
            vertex_indices: dict[tuple[int, ...], int] = {}
            minimum_bmac = stress_bmac
            case_status = "SUCCESS"
            for subset in mobius_vertices(coalition):
                tracked = track_action(gamma, empty_modes, empty_index, subset)
                if tracked is None:
                    case_status = "MODE_TRACK_FAIL"
                    break
                target_index, action_bmac, homotopy_nodes, status = tracked
                if status != "SUCCESS" or target_index < 0:
                    case_status = status
                    break
                point, target_modes = get(gamma, subset, 1.0 if subset else 0.0)
                if target_modes is None:
                    case_status = "EIG_FAIL"
                    break
                value = complex(target_modes.eigenvalues[target_index])
                values[subset] = value
                damping[subset] = damping_ratio(value)
                modes_for_vertex[subset] = target_modes
                vertex_indices[subset] = target_index
                minimum_bmac = min(minimum_bmac, action_bmac)
                track_rows.append(
                    {
                        **base,
                        "subset": subset_key(subset),
                        "lambda_real": float(value.real),
                        "lambda_imag": float(value.imag),
                        "frequency_hz": float(value.imag / (2.0 * math.pi)),
                        "damping_ratio": damping[subset],
                        "BMAC": float(min(stress_bmac, action_bmac)),
                        "stress_refinement_depth": stress_depth,
                        "action_homotopy_nodes": homotopy_nodes,
                        "kappa_lambda": float(target_modes.condition_numbers[target_index]),
                        "residual": float(target_modes.residuals[target_index]),
                        "nearest_modal_separation_per_s": float(target_modes.nearest_separations[target_index]),
                        "status": "SUCCESS",
                    }
                )
            if case_status != "SUCCESS":
                externality_rows.append({**base, "status": case_status, "minimum_BMAC": minimum_bmac})
                continue
            delta = complex_externality(values, coalition)
            delta_zeta = coalition_externality(damping, coalition)
            lower = stable_lower_order_prediction(values, coalition)
            full = values[tuple(sorted(coalition))]
            closure = full - (lower + delta)
            lambda0 = values[()]
            rho0 = delta.real / (-lambda0.real) if lambda0.real < 0.0 else math.nan
            rho_comp = delta.real / (-lower.real) if lower.real < 0.0 else math.nan
            lower_zeta = damping_ratio(lower)
            full_zeta = damping_ratio(full)
            hard_failure = bool(lower.real < 0.0 and full.real >= 0.0)
            screen_failure = bool(lower_zeta >= screen and full_zeta < screen)
            pole_values = {subset: get_pole(gamma, subset) for subset in mobius_vertices(coalition)}
            delta_phi_poles = coalition_externality(pole_values, coalition)
            analysis_status = "LOWER_ORDER_UNSTABLE" if lower.real >= 0.0 else "SUCCESS"
            externality_rows.append(
                {
                    **base,
                    "status": analysis_status,
                    "minimum_BMAC": float(minimum_bmac),
                    "lambda_empty_real": float(lambda0.real),
                    "lambda_empty_imag": float(lambda0.imag),
                    "lambda_full_real": float(full.real),
                    "lambda_full_imag": float(full.imag),
                    "Delta_lambda_real": float(delta.real),
                    "Delta_lambda_imag": float(delta.imag),
                    "Delta_lambda_absolute": float(abs(delta)),
                    "Delta_zeta": float(delta_zeta),
                    "lambda_lower_order_real": float(lower.real),
                    "lambda_lower_order_imag": float(lower.imag),
                    "zeta_lower_order": float(lower_zeta),
                    "zeta_full": float(full_zeta),
                    "rho0": float(rho0),
                    "rho_comp": float(rho_comp),
                    "composition_closure_real": float(closure.real),
                    "composition_closure_imag": float(closure.imag),
                    "hard_composition_failure": hard_failure,
                    "screen_composition_failure": screen_failure,
                    "stress_conditioned_material_delta_zeta": bool(abs(delta_zeta) >= 0.0025),
                    "destabilizing_material_delta_zeta": bool(delta_zeta <= -0.0025 and delta.real > 0.0),
                    "DeltaPhi_poles": float(delta_phi_poles),
                    "empty_kappa_lambda": float(empty_modes.condition_numbers[empty_index]),
                    "full_kappa_lambda": float(modes_for_vertex[tuple(sorted(coalition))].condition_numbers[vertex_indices[tuple(sorted(coalition))]]),
                    "empty_nearest_modal_separation_per_s": float(empty_modes.nearest_separations[empty_index]),
                }
            )
    return {
        "seed_id": seed["seed_id"],
        "tracks": track_rows,
        "externalities": externality_rows,
        "builds": build_rows,
        "wall_time_s": time.perf_counter() - started,
        "physical_build_count": len(cache),
        "pole_factor_count": len(pole_cache),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    limits_path = ARTIFACT / "baseline_calibration" / "E05C_BASELINE_STRESS_LIMITS.parquet"
    if not limits_path.is_file():
        raise FileNotFoundError("baseline-only stress limits must be frozen before coalition campaign")
    limit_manifest = PREREG / "E05C_STRESS_ENDPOINT_FREEZE.json"
    if not limit_manifest.is_file():
        raise RuntimeError("E05C_STRESS_ENDPOINT_FREEZE.json is required before coalition outcomes")
    candidates = load_json(PREREG / "E05C_CANDIDATE_FREEZE.json")["candidates"]
    mode_freeze = load_json(PREREG / "E05C_MODE_FAMILY_FREEZE.json")
    locks = {record["candidate_id"]: record for record in mode_freeze["mode_families"]}
    seeds = load_json(PREREG / "E05C_BASE_HOLDOUT_FREEZE.json")["operating_points"]
    limits = pd.read_parquet(limits_path).set_index("seed_id")
    stress = load_json(PREREG / "E05C_STRESS_RULE_FREEZE.json")
    numerical = load_json(PREREG / "E05C_NUMERICAL_FREEZE.json")
    gates = load_json(PREREG / "E05C_GATE_FREEZE.json")
    development = reference_modes()
    tasks = [
        {
            "seed": seed,
            "gamma_end": float(limits.loc[seed["seed_id"], "gamma_end"]),
            "candidates": candidates,
            "locks": locks,
            "stress": stress,
            "numerical": numerical,
            "gates": gates,
            "development_modes": development,
        }
        for seed in seeds
    ]
    started = time.perf_counter()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(worker, task): task["seed"]["seed_id"] for task in tasks}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(f"{count}/8 {result['seed_id']} {result['physical_build_count']} builds {result['wall_time_s']:.1f}s", flush=True)
    results.sort(key=lambda value: value["seed_id"])
    tracks = pd.DataFrame([row for result in results for row in result["tracks"]])
    externalities = pd.DataFrame([row for result in results for row in result["externalities"]])
    builds = pd.DataFrame([row for result in results for row in result["builds"]])
    output = ARTIFACT / "holdout"
    tracks.to_parquet(output / "e05c_mode_tracks.parquet", index=False)
    externalities.to_parquet(output / "e05c_connected_pole_externalities.parquet", index=False)
    externalities.to_parquet(output / "e05c_margin_consumption.parquet", index=False)
    externalities[
        [
            "seed_id", "tau", "gamma", "candidate_id", "coalition_key", "order", "modal_family_id",
            "lambda_lower_order_real", "lambda_lower_order_imag", "lambda_full_real", "lambda_full_imag",
            "zeta_lower_order", "zeta_full", "hard_composition_failure", "screen_composition_failure", "status",
        ]
    ].to_parquet(output / "e05c_composition_failures.parquet", index=False)
    tracks[
        [
            "seed_id", "tau", "gamma", "candidate_id", "coalition_key", "modal_family_id", "subset",
            "lambda_real", "lambda_imag", "BMAC", "kappa_lambda", "residual", "nearest_modal_separation_per_s", "status",
        ]
    ].to_parquet(output / "e05c_conditioning.parquet", index=False)
    externalities[
        ["seed_id", "tau", "gamma", "candidate_id", "coalition_key", "order", "DeltaPhi_poles", "status"]
    ].to_parquet(output / "e05c_dynamic_pole_externality.parquet", index=False)
    builds.to_parquet(output / "e05c_failure_ledger.parquet", index=False)
    compute = {
        "seed_count": 8,
        "tau_count": 9,
        "candidate_count": 15,
        "logical_candidate_vertex_requests": 8 * 9 * (10 * 4 + 5 * 8),
        "unique_endpoint_subset_count": len(union_vertices([parse_subset_key(record["coalition_key"]) for record in candidates])),
        "physical_build_count": sum(result["physical_build_count"] for result in results),
        "eigensolve_count": int((builds["maximum_eigenpair_residual"].notna()).sum()),
        "pole_factor_count": sum(result["pole_factor_count"] for result in results),
        "successful_externality_rows": int((externalities["status"] == "SUCCESS").sum()),
        "externality_rows": len(externalities),
        "mode_track_rows": len(tracks),
        "build_status_counts": builds["status"].value_counts(dropna=False).to_dict(),
        "workers": args.workers,
        "wall_time_s": time.perf_counter() - started,
    }
    (ARTIFACT / "logs" / "e05c_campaign_compute.json").write_text(
        json.dumps(compute, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(compute, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
