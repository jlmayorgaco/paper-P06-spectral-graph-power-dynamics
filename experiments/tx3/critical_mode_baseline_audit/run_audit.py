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

import andes
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E03), str(E05D)]

from campaign_core import action_vector, coalition_externality  # noqa: E402
from weak_grid import (  # noqa: E402
    CASE,
    GeneralizedModes,
    complex_externality,
    critical_mode_index,
    damping_ratio,
    eigenvector_hash,
    evaluate_weak_grid_operator,
    generalized_modes,
    mobius_vertices,
    parse_subset_key,
    stable_lower_order_prediction,
    subset_key,
    track_one_hungarian,
    union_vertices,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "critical_mode_baseline_audit"
E05D_PREREG = ROOT / "artifacts" / "tx3" / "E05D_weak_grid" / "preregistration"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def state_names() -> list[str]:
    system = andes.load(str(CASE), setup=True, no_output=True)
    if not system.PFlow.run() or not system.PFlow.converged:
        raise RuntimeError("state-name reference power flow failed")
    system.TDS.init()
    return [str(value) for value in system.dae.x_name]


def state_group(name: str) -> str:
    parts = name.split()
    state = parts[0] if parts else ""
    model = parts[1] if len(parts) > 1 else ""
    if model == "GENROU" and state in {"delta", "omega"}:
        return "synchronous_rotor"
    if model == "GENROU":
        return "synchronous_electromagnetic"
    if model == "GAST":
        return "turbine_governor"
    if model == "SEXS":
        return "excitation"
    if model in {"PLL2", "BusFreq"}:
        return "pll"
    if model == "REGCP1":
        return "converter_current_command"
    if model == "REECB1":
        return "converter_electrical_control"
    if model == "REPCA1":
        return "plant_control"
    return "other"


def mode_anatomy(modes: GeneralizedModes, index: int, names: list[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    left = modes.left[:, index]
    right = modes.right[:, index]
    raw_participation = np.abs(left.conj() * right)
    participation = raw_participation / max(float(raw_participation.sum()), 1e-30)
    denominator = np.vdot(left, modes.mass * right)
    raw_endogeneity = left.conj() * modes.mass * right / denominator
    endogeneity_abs = np.abs(raw_endogeneity)
    endogeneity_abs /= max(float(endogeneity_abs.sum()), 1e-30)
    rows = []
    for state_index, name in enumerate(names):
        rows.append({
            "state_index": state_index, "state_name": name, "state_group": state_group(name),
            "participation_abs_normalized": float(participation[state_index]),
            "endogeneity_real": float(raw_endogeneity[state_index].real),
            "endogeneity_imag": float(raw_endogeneity[state_index].imag),
            "endogeneity_abs_normalized": float(endogeneity_abs[state_index]),
        })
    fine = pd.DataFrame(rows).groupby("state_group")[["participation_abs_normalized", "endogeneity_abs_normalized"]].sum()
    def score(groups: set[str], column: str) -> float:
        return float(fine.reindex(sorted(groups), fill_value=0.0)[column].sum())
    synchronous_groups = {"synchronous_rotor", "synchronous_electromagnetic", "turbine_governor", "excitation"}
    current_groups = {"converter_current_command", "converter_electrical_control", "plant_control"}
    broad = {
        "synchronous_electromechanical": score(synchronous_groups, "participation_abs_normalized"),
        "PLL_dominated": score({"pll"}, "participation_abs_normalized"),
        "current_control_dominated": score(current_groups, "participation_abs_normalized"),
        "other": score({"other"}, "participation_abs_normalized"),
    }
    endogeneity = {
        "synchronous_electromechanical": score(synchronous_groups, "endogeneity_abs_normalized"),
        "PLL_dominated": score({"pll"}, "endogeneity_abs_normalized"),
        "current_control_dominated": score(current_groups, "endogeneity_abs_normalized"),
        "other": score({"other"}, "endogeneity_abs_normalized"),
    }
    ordered = sorted(broad.items(), key=lambda item: item[1], reverse=True)
    if ordered[0][1] >= 0.60 and ordered[0][0] != "other":
        classification = ordered[0][0]
    elif len(ordered) > 1 and ordered[0][1] >= 0.20 and ordered[1][1] >= 0.20:
        classification = "hybrid"
    else:
        classification = "other"
    summary = {
        "classification": classification, "classification_threshold": 0.60,
        "participation_scores": broad, "endogeneity_scores": endogeneity,
        "top_states": sorted(rows, key=lambda row: row["participation_abs_normalized"], reverse=True)[:12],
    }
    return rows, summary


def worker(task: dict[str, Any]) -> dict[str, Any]:
    seed = task["seed"]
    candidates = task["candidates"]
    names = task["state_names"]
    bmac_min = float(task["protocol"]["BMAC_minimum"])
    residual_max = float(task["protocol"]["eigenpair_residual_maximum"])
    started = time.perf_counter()
    cache: dict[tuple[tuple[int, ...], float], tuple[Any, GeneralizedModes | None, str]] = {}
    build_rows: list[dict[str, Any]] = []
    requests = 0

    def key(subset: tuple[int, ...], beta: float) -> tuple[tuple[int, ...], float]:
        return ((), 0.0) if beta == 0.0 or not subset else (tuple(subset), round(float(beta), 12))

    def get(subset: tuple[int, ...] = (), beta: float = 0.0) -> tuple[Any, GeneralizedModes | None, str]:
        nonlocal requests
        requests += 1
        cache_key = key(subset, beta)
        if cache_key in cache:
            return cache[cache_key]
        point = evaluate_weak_grid_operator(seed, 1.0, action_vector({f"A{value}": float(beta) for value in subset}))
        modes = None
        status = point.status
        eig_exception = None
        if status == "SUCCESS":
            try:
                modes = generalized_modes(point)
                if float(np.max(modes.residuals)) > residual_max:
                    status = "EIGENPAIR_RESIDUAL_FAIL"
            except Exception as exc:
                status = "EIG_FAIL"; eig_exception = repr(exc)
        cache[cache_key] = (point, modes, status)
        build_rows.append({
            "seed_id": seed["seed_id"], "kappa_grid": 1.0, "subset": subset_key(subset),
            "beta": float(beta), "status": status, "eig_exception": eig_exception,
            "maximum_eigenpair_residual": float(np.max(modes.residuals)) if modes is not None else None,
            **point.metadata,
        })
        return cache[cache_key]

    empty_point, empty_modes, empty_status = get()
    if empty_status != "SUCCESS" or empty_modes is None:
        raise RuntimeError(f"{seed['seed_id']} empty baseline failed: {empty_status}")
    empty_index = critical_mode_index(empty_modes)
    empty_value = complex(empty_modes.eigenvalues[empty_index])
    anatomy_rows, anatomy = mode_anatomy(empty_modes, empty_index, names)
    for row in anatomy_rows:
        row["seed_id"] = seed["seed_id"]
    mode_row = {
        "seed_id": seed["seed_id"], "kappa_grid": 1.0, "mode_index": empty_index,
        "lambda_real": float(empty_value.real), "lambda_imag": float(empty_value.imag),
        "frequency_hz": float(empty_value.imag / (2.0 * math.pi)),
        "damping_ratio": damping_ratio(empty_value), "BMAC": 1.0,
        "eigenpair_residual": float(empty_modes.residuals[empty_index]),
        "kappa_lambda": float(empty_modes.condition_numbers[empty_index]),
        "nearest_modal_separation_per_s": float(empty_modes.nearest_separations[empty_index]),
        "left_right_eigenvectors_sha256": eigenvector_hash(empty_modes, empty_index),
        "classification": anatomy["classification"],
        "synchronous_participation": anatomy["participation_scores"]["synchronous_electromechanical"],
        "pll_participation": anatomy["participation_scores"]["PLL_dominated"],
        "current_control_participation": anatomy["participation_scores"]["current_control_dominated"],
        "other_participation": anatomy["participation_scores"]["other"],
        "synchronous_endogeneity": anatomy["endogeneity_scores"]["synchronous_electromechanical"],
        "pll_endogeneity": anatomy["endogeneity_scores"]["PLL_dominated"],
        "current_control_endogeneity": anatomy["endogeneity_scores"]["current_control_dominated"],
        "other_endogeneity": anatomy["endogeneity_scores"]["other"],
        "top_states_json": json.dumps(anatomy["top_states"], sort_keys=True),
        "SCR36": float(empty_point.metadata["SCR36"]), "SCR37": float(empty_point.metadata["SCR37"]),
        "SCR38": float(empty_point.metadata["SCR38"]), "min_SCR": float(empty_point.metadata["min_SCR"]),
        "status": "SUCCESS",
    }

    def track(subset: tuple[int, ...]) -> tuple[int, float, int, str]:
        if not subset:
            return empty_index, 1.0, 1, "SUCCESS"
        last_status = "MODE_TRACK_FAIL"
        for count in (5, 9, 17):
            previous_modes = empty_modes; previous_index = empty_index; minimum = 1.0; success = True
            for beta in np.linspace(0.0, 1.0, count)[1:]:
                _, modes, status = get(subset, float(beta))
                if status != "SUCCESS" or modes is None:
                    success = False; last_status = "NONSMOOTH_REGIME" if status == "LIMITER_ACTIVE" else status; break
                match = track_one_hungarian(previous_modes, modes, previous_index)
                if match is None or match[1] < bmac_min:
                    success = False; last_status = "MODE_TRACK_FAIL"; break
                previous_index, score = match; previous_modes = modes; minimum = min(minimum, score)
            if success:
                return previous_index, minimum, count, "SUCCESS"
        return -1, 0.0, 17, last_status

    coalitions = [parse_subset_key(record["coalition_key"]) for record in candidates]
    tracked: dict[tuple[int, ...], dict[str, Any]] = {}
    subset_rows: list[dict[str, Any]] = []
    for subset in union_vertices(coalitions):
        target_index, bmac, nodes, status = track(subset)
        if status != "SUCCESS":
            tracked[subset] = {"status": status}
            subset_rows.append({"seed_id": seed["seed_id"], "subset": subset_key(subset), "status": status})
            continue
        _, modes, endpoint_status = get(subset, 1.0 if subset else 0.0)
        if endpoint_status != "SUCCESS" or modes is None:
            tracked[subset] = {"status": endpoint_status}; continue
        value = complex(modes.eigenvalues[target_index])
        row = {
            "seed_id": seed["seed_id"], "kappa_grid": 1.0, "subset": subset_key(subset),
            "order": len(subset), "mode_index": int(target_index),
            "lambda_real": float(value.real), "lambda_imag": float(value.imag),
            "frequency_hz": float(value.imag / (2.0 * math.pi)), "damping_ratio": damping_ratio(value),
            "BMAC": float(bmac), "action_homotopy_nodes": int(nodes),
            "eigenpair_residual": float(modes.residuals[target_index]),
            "kappa_lambda": float(modes.condition_numbers[target_index]),
            "nearest_modal_separation_per_s": float(modes.nearest_separations[target_index]),
            "status": "SUCCESS",
        }
        tracked[subset] = row; subset_rows.append(row)

    externality_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        coalition = parse_subset_key(candidate["coalition_key"])
        required = mobius_vertices(coalition)
        base = {
            "seed_id": seed["seed_id"], "kappa_grid": 1.0,
            "candidate_id": candidate["candidate_id"], "coalition": candidate["coalition_key"],
            "order": int(candidate["order"]), "critical_mode_classification": anatomy["classification"],
        }
        failure = next((tracked[subset]["status"] for subset in required if tracked[subset]["status"] != "SUCCESS"), None)
        if failure is not None:
            externality_rows.append({**base, "status": failure}); continue
        values = {subset: complex(tracked[subset]["lambda_real"], tracked[subset]["lambda_imag"]) for subset in required}
        damping = {subset: float(tracked[subset]["damping_ratio"]) for subset in required}
        delta = complex_externality(values, coalition)
        delta_zeta = coalition_externality(damping, coalition)
        lower = stable_lower_order_prediction(values, coalition)
        full = values[tuple(sorted(coalition))]
        closure = full - (lower + delta)
        rho_comp = delta.real / (-lower.real) if lower.real < 0.0 else math.nan
        rho_baseline = delta.real / (-empty_value.real) if empty_value.real < 0.0 else math.nan
        externality_rows.append({
            **base, "lambda_empty_real": float(empty_value.real), "lambda_empty_imag": float(empty_value.imag),
            "lambda_full_real": float(full.real), "lambda_full_imag": float(full.imag),
            "Delta_lambda_real": float(delta.real), "Delta_lambda_imag": float(delta.imag),
            "Delta_lambda_absolute": float(abs(delta)), "Delta_zeta": float(delta_zeta),
            "lambda_lower_real": float(lower.real), "lambda_lower_imag": float(lower.imag),
            "zeta_lower": damping_ratio(lower), "zeta_full": damping_ratio(full),
            "rho_comp": float(rho_comp), "rho_baseline": float(rho_baseline),
            "hard_failure": bool(lower.real < 0.0 and full.real >= 0.0),
            "screen_failure": bool(damping_ratio(lower) >= 0.05 and damping_ratio(full) < 0.05),
            "minimum_BMAC": float(min(tracked[subset]["BMAC"] for subset in required)),
            "composition_closure_real": float(closure.real), "composition_closure_imag": float(closure.imag),
            "status": "LOWER_ORDER_UNSTABLE" if lower.real >= 0.0 else "SUCCESS",
        })
    return {
        "seed_id": seed["seed_id"], "mode": mode_row, "anatomy": anatomy_rows,
        "subsets": subset_rows, "externalities": externality_rows, "builds": build_rows,
        "left": empty_modes.left[:, empty_index], "right": empty_modes.right[:, empty_index],
        "physical_build_count": len(cache), "cache_requests": requests,
        "wall_time_s": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    protocol_path = ARTIFACT / "protocol" / "CRITICAL_MODE_BASELINE_AUDIT_FREEZE.json"
    if not protocol_path.is_file():
        raise FileNotFoundError("audit protocol freeze is required")
    if (ARTIFACT / "tables" / "critical_mode_externalities.parquet").exists():
        raise RuntimeError("audit results already exist")
    protocol = load_json(protocol_path)
    seeds = load_json(E05D_PREREG / "E05D_HOLDOUT_FREEZE.json")["operating_points"]
    candidates = load_json(E05D_PREREG / "E05D_CANDIDATE_FREEZE.json")["candidates"]
    names = state_names()
    tasks = [{"seed": seed, "candidates": candidates, "state_names": names, "protocol": protocol} for seed in seeds]
    started = time.perf_counter(); results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(worker, task): task["seed"]["seed_id"] for task in tasks}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result(); results.append(result)
            print(f"{count}/8 {result['seed_id']} {result['physical_build_count']} builds {result['wall_time_s']:.1f}s", flush=True)
    results.sort(key=lambda result: result["seed_id"])
    tables = ARTIFACT / "tables"; logs = ARTIFACT / "logs"; tables.mkdir(parents=True, exist_ok=True); logs.mkdir(parents=True, exist_ok=True)
    modes = pd.DataFrame([result["mode"] for result in results])
    anatomy = pd.DataFrame([row for result in results for row in result["anatomy"]])
    subsets = pd.DataFrame([row for result in results for row in result["subsets"]])
    externalities = pd.DataFrame([row for result in results for row in result["externalities"]])
    builds = pd.DataFrame([row for result in results for row in result["builds"]])
    modes.to_parquet(tables / "critical_modes.parquet", index=False)
    anatomy.to_parquet(tables / "critical_mode_state_anatomy.parquet", index=False)
    subsets.to_parquet(tables / "critical_mode_subset_poles.parquet", index=False)
    externalities.to_parquet(tables / "critical_mode_externalities.parquet", index=False)
    builds.to_parquet(tables / "audit_failure_ledger.parquet", index=False)
    np.savez_compressed(ARTIFACT / "critical_mode_eigenvectors.npz", **{
        f"{result['seed_id']}_{side}": result[side] for result in results for side in ("left", "right")
    })
    requests = sum(result["cache_requests"] for result in results); builds_count = sum(result["physical_build_count"] for result in results)
    summary = {
        "scope": "descriptive post-mortem at kappa_grid=1; no claim gate",
        "seed_count": 8, "candidate_count": 15, "subset_count_per_seed": len(subsets) // 8,
        "physical_build_count": builds_count, "eigensolve_count": int(builds["maximum_eigenpair_residual"].notna().sum()),
        "cache_requests": requests, "cache_hits": requests - builds_count,
        "cache_hit_rate": (requests - builds_count) / requests,
        "mode_track_success_rate": float(subsets["status"].eq("SUCCESS").mean()),
        "externality_success_rate": float(externalities["status"].eq("SUCCESS").mean()),
        "build_status_counts": builds["status"].value_counts().to_dict(),
        "wall_time_s": time.perf_counter() - started,
    }
    (logs / "audit_compute.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
