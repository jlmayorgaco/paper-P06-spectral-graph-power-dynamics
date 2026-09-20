"""Preregistered exact reduced-DAE TX4 robustness campaign.

The script evaluates every one of the 16 target-family portfolios for every
condition in the deterministic, QMC, and Monte-Carlo designs.  Conditions are
independent and resume-safe at the row level: a completed CSV is never
silently replaced by a smaller sample.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research" / "ias2026_last_validation" / "code" / "python"
sys.path.insert(0, str(SOURCE))

from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import Ieee39Network, load_network  # noqa: E402


TARGET_BUSES = (30, 33, 35, 37)
PORTFOLIOS = [
    tuple(TARGET_BUSES[i] for i in range(4) if mask & (1 << i))
    for mask in range(16)
]
PORTFOLIO_LABELS = {p: ("BASE" if not p else "+".join(map(str, p))) for p in PORTFOLIOS}
PARAM_NAMES = ("g", "k", "t", "h", "epsilon", "damping", "inertia")
BOUNDS = {
    "g": (0.020, 0.250),
    "k": (0.750, 1.750),
    "t": (0.750, 1.750),
    "h": (0.500, 1.500),
    "epsilon": (0.900, 1.100),
    "damping": (0.000, 0.200),
    "inertia": (0.800, 1.200),
}
NOMINAL = {"g": 0.03625, "k": 1.425, "t": 1.5, "h": 1.0, "epsilon": 1.0, "damping": 0.0, "inertia": 1.0}
THRESHOLD = 1e-8
GAUGE_TOL = 1e-4
MODE_LOW_HZ = 0.3
MODE_HIGH_HZ = 1.5

FIELDS = [
    "campaign", "condition_id", "seed", "portfolio", "replaced_count",
    *PARAM_NAMES, "status", "error_type", "error_message", "states",
    "norm_f", "norm_g", "gz_condition", "ell_H4", "q_H4",
    "alpha_all", "alpha_EM", "critical_frequency_hz", "stable_all",
    "stable_EM", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE", "delta_H4",
    "eta_H4", "H4_alpha_all", "proper_max_alpha_all", "proper_stable_count",
    "proper_count", "condition_runtime_s",
]


def network_for_epsilon(epsilon: float) -> Ieee39Network:
    base = load_network()
    loads = {bus: value * float(epsilon) for bus, value in base.loads.items()}
    return replace(base, loads=loads)


def _float(value: float | int | None) -> float | None:
    return None if value is None else float(value)


def _one_portfolio(params: dict[str, float], portfolio: tuple[int, ...], campaign: str, condition_id: str, seed: int) -> dict[str, Any]:
    start = time.perf_counter()
    base = {
        "campaign": campaign,
        "condition_id": condition_id,
        "seed": seed,
        "portfolio": PORTFOLIO_LABELS[portfolio],
        "replaced_count": len(portfolio),
        **{key: float(params[key]) for key in PARAM_NAMES},
        "status": "OK",
        "error_type": "",
        "error_message": "",
        "states": None,
        "norm_f": None,
        "norm_g": None,
        "gz_condition": None,
        "ell_H4": None,
        "q_H4": None,
        "alpha_all": None,
        "alpha_EM": None,
        "critical_frequency_hz": None,
        "stable_all": None,
        "stable_EM": None,
        "H4_PRESENT": None,
        "EXACT_H4": None,
        "NONCOMPOSABLE": False,
        "delta_H4": None,
        "eta_H4": None,
        "H4_alpha_all": None,
        "proper_max_alpha_all": None,
        "proper_stable_count": None,
        "proper_count": 15,
        "condition_runtime_s": None,
    }
    try:
        network = network_for_epsilon(params["epsilon"])
        converter = ConverterParameters(
            voltage_control=True,
            voltage_gain=params["g"],
            voltage_leak=0.05,
        )
        converters = None
        if 37 in portfolio:
            converters = {
                37: ConverterParameters(
                    voltage_control=True,
                    voltage_gain=params["g"] * params["h"],
                    voltage_leak=0.05,
                )
            }
        plan = ReplacementPlan.of(
            {bus: 1.0 for bus in portfolio},
            q_policy="matched",
            device="gfl",
            machine_damping=params["damping"],
        )
        case = solve_case(
            plan,
            network=network,
            converter=converter,
            converters=converters,
            machine_scaling={"ka": params["k"], "ta": params["t"]},
            machine_services={"inertia": params["inertia"]},
            tol=1e-9,
        )
        jac = central_difference_jacobians(
            case.dae, case.equilibrium.x, case.equilibrium.z, case.equilibrium.theta
        )
        eig = np.linalg.eigvals(case.system.A)
        transverse = eig[np.abs(eig) >= GAUGE_TOL]
        if transverse.size == 0:
            transverse = eig
        frequency = np.abs(transverse.imag) / (2.0 * np.pi)
        em = transverse[(frequency >= MODE_LOW_HZ) & (frequency <= MODE_HIGH_HZ)]
        alpha_all = float(np.max(eig.real))
        alpha_em = float(np.max(em.real)) if em.size else float("nan")
        critical = transverse[int(np.argmax(transverse.real))]
        # ell_H4 is an algebraic local-regularity diagnostic; q_H4 is a
        # modal distance-to-marginality diagnostic. Neither is a physical radius.
        sigma = np.linalg.svd(jac.gz, compute_uv=False)
        ell = float(np.min(sigma))
        q = float(np.min(np.abs(transverse.real)))
        base.update(
            {
                "states": int(case.n_states),
                "norm_f": float(case.equilibrium.norm_f),
                "norm_g": float(case.equilibrium.norm_g),
                "gz_condition": float(case.gz_condition),
                "ell_H4": ell,
                "q_H4": q,
                "alpha_all": alpha_all,
                "alpha_EM": alpha_em,
                "critical_frequency_hz": float(abs(critical.imag) / (2.0 * np.pi)),
                "stable_all": bool(alpha_all <= THRESHOLD),
                "stable_EM": bool(np.isfinite(alpha_em) and alpha_em <= THRESHOLD),
            }
        )
    except Exception as exc:  # retain failed conditions as NONCOMPOSABLE
        base.update(
            {
                "status": "NONCOMPOSABLE",
                "error_type": type(exc).__name__,
                "error_message": str(exc)[:500],
                "NONCOMPOSABLE": True,
            }
        )
    base["condition_runtime_s"] = time.perf_counter() - start
    return base


def attach_condition_summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    h4 = next((row for row in rows if row["portfolio"] == "30+33+35+37"), None)
    proper = [row for row in rows if row["portfolio"] != "30+33+35+37"]
    good_proper = [row for row in proper if row["status"] == "OK" and row["alpha_all"] is not None]
    h4_present = bool(h4 and h4["status"] == "OK" and h4["alpha_all"] is not None and h4["alpha_all"] > THRESHOLD)
    proper_max = max((float(row["alpha_all"]) for row in good_proper), default=float("nan"))
    proper_stable = sum(bool(row["stable_all"]) for row in good_proper)
    exact = bool(h4_present and len(good_proper) == 15 and proper_stable == 15)
    h4_alpha = None if not h4 or h4["alpha_all"] is None else float(h4["alpha_all"])
    delta = None if h4_alpha is None or not np.isfinite(proper_max) else h4_alpha - proper_max
    eta = None if not np.isfinite(proper_max) else -proper_max
    for row in rows:
        row.update(
            {
                "H4_PRESENT": h4_present,
                "EXACT_H4": exact,
                "H4_alpha_all": h4_alpha,
                "proper_max_alpha_all": None if not np.isfinite(proper_max) else proper_max,
                "proper_stable_count": proper_stable,
                "delta_H4": delta,
                "eta_H4": eta,
            }
        )
    return rows


def evaluate_condition(item: tuple[str, str, int, dict[str, float]]) -> list[dict[str, Any]]:
    campaign, condition_id, seed, params = item
    rows = [_one_portfolio(params, portfolio, campaign, condition_id, seed) for portfolio in PORTFOLIOS]
    return attach_condition_summaries(rows)


def nominal() -> dict[str, float]:
    return dict(NOMINAL)


def make_conditions() -> list[tuple[str, str, int, dict[str, float]]]:
    conditions: list[tuple[str, str, int, dict[str, float]]] = []
    conditions.append(("NOMINAL", "nominal", 0, nominal()))
    for name in PARAM_NAMES:
        lo, hi = BOUNDS[name]
        for index, value in enumerate(np.linspace(lo, hi, 201)):
            p = nominal()
            p[name] = float(value)
            conditions.append(("1D", f"1d_{name}_{index:03d}", 0, p))
    pairs = [("g", "k"), ("g", "t"), ("g", "h"), ("g", "epsilon"), ("k", "epsilon"), ("t", "epsilon")]
    for first, second in pairs:
        lo1, hi1 = BOUNDS[first]
        lo2, hi2 = BOUNDS[second]
        values1 = np.linspace(lo1, hi1, 41)
        values2 = np.linspace(lo2, hi2, 41)
        for i, a in enumerate(values1):
            for j, b in enumerate(values2):
                p = nominal()
                p[first] = float(a)
                p[second] = float(b)
                conditions.append(("2D", f"2d_{first}x{second}_{i:02d}_{j:02d}", 0, p))
    sob = qmc.Sobol(d=len(PARAM_NAMES), scramble=True, seed=20260920)
    qmc_points = sob.random_base2(m=12)
    for index, point in enumerate(qmc_points):
        p = {name: BOUNDS[name][0] + float(point[i]) * (BOUNDS[name][1] - BOUNDS[name][0]) for i, name in enumerate(PARAM_NAMES)}
        conditions.append(("QMC", f"qmc_{index:04d}", 20260920, p))
    rng = np.random.default_rng(20260921)
    mc_points = rng.random((5000, len(PARAM_NAMES)))
    for index, point in enumerate(mc_points):
        p = {name: BOUNDS[name][0] + float(point[i]) * (BOUNDS[name][1] - BOUNDS[name][0]) for i, name in enumerate(PARAM_NAMES)}
        conditions.append(("MC", f"mc_{index:04d}", 20260921, p))
    return conditions


def write_rows(path: Path, rows: list[dict[str, Any]], mode: str = "a") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and mode == "a"
    with path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, min(18, (os.cpu_count() or 2) - 2)))
    parser.add_argument("--limit", type=int, default=0, help="engineering-only prefix limit; never used for the preregistered run")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv")
    args = parser.parse_args()
    conditions = make_conditions()
    if args.limit:
        conditions = conditions[: args.limit]
    if args.output.exists():
        args.output.unlink()
    start = time.perf_counter()
    completed = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(evaluate_condition, item) for item in conditions]
        for future in as_completed(futures):
            rows = future.result()
            write_rows(args.output, rows)
            completed += 1
            if completed % 25 == 0 or completed == len(conditions):
                elapsed = time.perf_counter() - start
                print(f"TX4_ROBUSTNESS_PROGRESS conditions={completed}/{len(conditions)} rows={completed*16} elapsed_s={elapsed:.1f}", flush=True)
    metadata = {
        "conditions": len(conditions),
        "rows": len(conditions) * 16,
        "workers": args.workers,
        "bounds": BOUNDS,
        "nominal": NOMINAL,
        "threshold": THRESHOLD,
        "qmc_seed": 20260920,
        "mc_seed": 20260921,
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "runtime_s": time.perf_counter() - start,
        "status": "COMPLETE" if not args.limit else "ENGINEERING_PREFIX_ONLY",
    }
    (args.output.parent / "TX4_ROBUSTNESS_MASTER_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    run()
