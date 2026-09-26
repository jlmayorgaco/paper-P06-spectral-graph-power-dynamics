"""E23 - Track-A validation gates 2 and 3.

Gate 2, closure continuation
    Follow the gauge-invariant closure distance ``d = min |lambda_i(Q) + 1|``
    from the empty portfolio up through every proper subset to the full core, and
    then through each repair. If the descriptor is meaningful it must approach
    zero only for the complete portfolio and move away again under repair.

Gate 3, mitigation frontier
    A single repair point is not an engineering answer. For a ladder of margin
    targets, find the minimum converter retune, the minimum synchronous condenser
    rating, and the resulting trade of controller change against synchronous MVA
    against retained PV.

Status of the result
    NUMERICAL OBSERVATION at one operating point, rho = 1.

Usage
    python experiments/E23_closure_and_frontier.py
"""

from __future__ import annotations

import argparse
import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_admittance import build_action_space

EXPERIMENT = "E23_closure_and_frontier"
CORE = (30, 33, 35, 37)
TARGETS = (-0.02, -0.05, -0.10, -0.126)
TUNED = ("kp_p", "ki_p", "kp_q", "ki_q", "tau_p", "kp_i", "ki_i")
BASE_CONVERTER = ConverterParameters()


def main() -> int:
    parser = argparse.ArgumentParser(description="Track-A gates 2 and 3")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    created = max(
        (m for m in eigen_analysis(flagship.system.A).modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    ).value
    total_pv = flagship.replaced_mw
    network = base.dae.network

    def closure(case, members) -> float:
        try:
            space = build_action_space(base, case, members)
            return float(abs(space.split(created)["closest_to_minus_one"] + 1.0))
        except (np.linalg.LinAlgError, ValueError):
            return float("nan")

    def abscissa(case) -> float:
        return assess(
            eigen_analysis(case.system.A), case.system.labels
        ).spectral_abscissa

    # ---- gate 2: continuation over the subset lattice ---------------------
    lattice = []
    for size in range(1, len(CORE) + 1):
        for members in combinations(CORE, size):
            case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
            lattice.append(
                {
                    "stage": "subset",
                    "members": "+".join(map(str, members)),
                    "size": size,
                    "spectral_abscissa": abscissa(case),
                    "closure_distance": closure(case, members),
                }
            )

    # ---- gate 2 continued: closure under each repair ----------------------
    def with_converter(converter) -> tuple[float, float]:
        case = solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE}), converter=converter
        )
        return abscissa(case), closure(case, CORE)

    def with_condenser(rating: float) -> tuple[float, float]:
        plan = ReplacementPlan.of(
            {b: 1.0 for b in CORE}, condenser={b: rating for b in CORE}
        )
        case = solve_case(plan)
        return abscissa(case), closure(case, CORE)

    repairs = []
    alpha, distance = with_converter(ConverterParameters(voltage_control=True))
    repairs.append(
        {
            "stage": "repair",
            "members": "voltage control",
            "size": 4,
            "spectral_abscissa": alpha,
            "closure_distance": distance,
        }
    )
    for scale in (0.75, 0.50, 0.43, 0.25):
        alpha, distance = with_converter(
            ConverterParameters(
                kp_p=BASE_CONVERTER.kp_p * scale, ki_p=BASE_CONVERTER.ki_p * scale
            )
        )
        repairs.append(
            {
                "stage": "repair",
                "members": f"P loop x{scale:g}",
                "size": 4,
                "spectral_abscissa": alpha,
                "closure_distance": distance,
            }
        )
    for rating in (0.05, 0.10, 0.25, 0.50):
        try:
            alpha, distance = with_condenser(rating)
        except Exception:
            continue
        repairs.append(
            {
                "stage": "repair",
                "members": f"condenser {rating:g}",
                "size": 4,
                "spectral_abscissa": alpha,
                "closure_distance": distance,
            }
        )
    continuation = pd.DataFrame(lattice + repairs)

    # ---- gate 3: mitigation frontier --------------------------------------
    baseline = np.array([getattr(BASE_CONVERTER, n) for n in TUNED])

    def build(vector: np.ndarray) -> ConverterParameters:
        return ConverterParameters(
            **dict(zip(TUNED, baseline * np.exp(vector), strict=True))
        )

    def retune_abscissa(vector: np.ndarray) -> float:
        return abscissa(
            solve_case(
                ReplacementPlan.of({b: 1.0 for b in CORE}), converter=build(vector)
            )
        )

    condenser_scan = []
    for rating in np.round(np.arange(0.05, 0.85, 0.05), 3):
        try:
            condenser_scan.append((float(rating), with_condenser(float(rating))[0]))
        except Exception:
            condenser_scan.append((float(rating), float("nan")))
    feasible_ratings = [r for r, a in condenser_scan if np.isfinite(a)]

    frontier = []
    for target in TARGETS:
        row = {"target_abscissa": target}
        # minimum-norm converter retune
        solution = minimize(
            lambda v: float(v @ v),
            np.full(len(TUNED), -0.2),
            method="SLSQP",
            bounds=[(-np.log(30.0), np.log(30.0))] * len(TUNED),
            constraints=[
                {"type": "ineq", "fun": lambda v, t=target: t - retune_abscissa(v)}
            ],
            options={"maxiter": 150, "ftol": 1e-9},
        )
        achieved = retune_abscissa(solution.x)
        row["retune_change_norm"] = (
            float(np.linalg.norm(solution.x)) if achieved <= target + 1e-6 else np.nan
        )
        row["retune_abscissa"] = achieved
        row["retune_pv_retained"] = total_pv
        # Minimum condenser rating. A bracketed root find is wrong here: small
        # ratings are INFEASIBLE, not merely insufficient, because the condenser
        # would have to carry the bus reactive output beyond its own rating. The
        # feasible set is scanned instead and the smallest working rating kept.
        working = [r for r, a in condenser_scan if np.isfinite(a) and a <= target]
        if working:
            rating = float(min(working))
            row["condenser_rating"] = rating
            row["condenser_mva"] = sum(rating * network.machines[b]["Sn"] for b in CORE)
        else:
            row["condenser_rating"] = np.nan
            row["condenser_mva"] = np.nan
        # cheapest single machine restoration that reaches the target
        best = None
        for bus in CORE:
            members = tuple(b for b in CORE if b != bus)
            case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
            if abscissa(case) <= target:
                mva = network.machines[bus]["Sn"]
                if best is None or mva < best[1]:
                    best = (bus, mva, case.replaced_mw)
        if best is not None:
            row["restore_bus"] = best[0]
            row["restore_mva"] = best[1]
            row["restore_pv_retained"] = best[2]
        frontier.append(row)
    frontier_table = pd.DataFrame(frontier)
    pd.DataFrame(condenser_scan, columns=["rating", "spectral_abscissa"]).to_csv(
        RESULTS / "tables" / f"{EXPERIMENT}_condenser_scan.csv", index=False
    )

    tables = RESULTS / "tables"
    continuation.to_csv(tables / f"{EXPERIMENT}_continuation.csv", index=False)
    frontier_table.to_csv(tables / f"{EXPERIMENT}_frontier.csv", index=False)

    subsets_only = continuation[continuation.stage == "subset"]
    proper = subsets_only[subsets_only["size"] < 4]
    full = subsets_only[subsets_only["size"] == 4].iloc[0]
    repairs_only = continuation[continuation.stage == "repair"]
    monotone = bool(
        full.closure_distance < proper.closure_distance.min() / 100.0
        and repairs_only[repairs_only.spectral_abscissa < 0].closure_distance.min()
        > full.closure_distance * 100.0
    )

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"core": list(CORE), "targets": list(TARGETS), "tuned": list(TUNED)},
    )
    manifest.finish(
        "SUCCESS" if monotone else "CLOSURE_NOT_DISCRIMINATING",
        closure_full_portfolio=float(full.closure_distance),
        closure_proper_subset_min=float(proper.closure_distance.min()),
        closure_repairs_min=float(
            repairs_only[repairs_only.spectral_abscissa < 0].closure_distance.min()
        ),
        frontier=frontier_table.to_dict("records"),
        feasible_condenser_ratings=feasible_ratings,
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("GATE 2  closure continuation, d = min |eig(Q) + 1| evaluated at the")
    print("        flagship crossing eigenvalue")
    print()
    print(continuation.to_string(index=False, float_format=lambda v: f"{v:12.6f}"))
    print()
    print(
        "  smallest closure distance over the proper subsets : %.4f"
        % proper.closure_distance.min()
    )
    print(
        "  closure distance of the full portfolio            : %.3e"
        % full.closure_distance
    )
    print(
        "  smallest closure distance over successful repairs : %.4f"
        % repairs_only[repairs_only.spectral_abscissa < 0].closure_distance.min()
    )
    print()
    print("GATE 3  mitigation frontier")
    print()
    print(frontier_table.to_string(index=False, float_format=lambda v: f"{v:12.4f}"))
    print()
    print(
        "GATE 2 PASSED, the closure descriptor is sharp and reversible"
        if monotone
        else "GATE 2 FAILED, the closure descriptor does not discriminate"
    )
    print(f"manifest -> {path}")
    return 0 if monotone else 1


if __name__ == "__main__":
    raise SystemExit(main())
