"""E34 - overnight 6. Mitigation Pareto frontier across four declared targets.

Not one target and not one cost. Four margins are declared before optimizing:

    -0.02, -0.05, -0.10 and the base case itself, alpha_base = -0.126478

and five strategy families are costed on axes that are reported SEPARATELY,
never collapsed into an arbitrary scalar:

    M1 converter control only
    M2 synchronous condenser only
    M3 converter control plus condenser
    M4 restore one synchronous machine
    M5 restore fractional synchronous capacity

Controller decision variables are physical, not raw gains: phase-locked-loop
natural frequency and damping, outer active and reactive bandwidths, and the
measurement-filter bandwidth.

The constrained quantity is the SPECTRAL ABSCISSA, which is what
``trackA_mitigation_definitions.yaml`` froze and what E21 used. The tracked
inter-area family envelope is reported beside it in every row, per the same
frozen rule that a mitigation number without both quantities is not a result.

The closure margin is measured before and after every repair, to ask whether a
successful repair moves the interaction spectrum away from -1 or merely moves an
eigenvalue.
"""

from __future__ import annotations

import time
from dataclasses import replace
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    SYSTEM_BASE_MVA,
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "E34_mitigation_frontier"
BAND_SAMPLES = 61
TARGETS = (-0.02, -0.05, -0.10, -0.126478)
TARGET_NAMES = ("-0.02", "-0.05", "-0.10", "alpha_base")
ZERO_MODE = 1e-3

#: Physical controller coordinates, with their nominal values.
PHYSICAL = ("pll_wn", "pll_zeta", "p_bandwidth", "q_bandwidth", "filter_bandwidth")
NOMINAL_PLL_WN = float(np.sqrt(ConverterParameters().ki_pll))
NOMINAL_PLL_ZETA = float(
    ConverterParameters().kp_pll / (2.0 * np.sqrt(ConverterParameters().ki_pll))
)


def converter_from_physical(vector):
    """Map log-deviations in physical coordinates to converter parameters."""

    base = ConverterParameters()
    wn = NOMINAL_PLL_WN * np.exp(vector[0])
    zeta = NOMINAL_PLL_ZETA * np.exp(vector[1])
    pll = ConverterParameters.pll_from_design(float(wn), float(zeta))
    p_scale = float(np.exp(vector[2]))
    q_scale = float(np.exp(vector[3]))
    filter_scale = float(np.exp(vector[4]))
    return replace(
        base,
        kp_pll=pll["kp_pll"],
        ki_pll=pll["ki_pll"],
        kp_p=base.kp_p * p_scale,
        ki_p=base.ki_p * p_scale,
        kp_q=base.kp_q * q_scale,
        ki_q=base.ki_q * q_scale,
        tau_p=base.tau_p / filter_scale,
    )


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="What does each mitigation cost, on separate axes, at four declared margins?",
        config={
            "targets": list(TARGETS),
            "constrained_quantity": "spectral abscissa (frozen definition)",
            "controller_coordinates": list(PHYSICAL),
            "cost_axes": [
                "pv_mw_forgone",
                "synchronous_mva_added",
                "controller_change_norm",
            ],
            "closure_metric": "frozen imaginary-axis m4 over the inter-area band",
        },
        workers=1,
    )
    started = time.time()
    network = load_network()
    base_case = solve_case(ReplacementPlan.of({}))
    flagship_plan = ReplacementPlan.of({b: 1.0 for b in CORE})
    flagship = solve_case(flagship_plan)
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base_case.system.A), base_case, nominal)
    total_pv_mw = float(flagship.replaced_mw)

    def measure(case):
        spectrum = eigen_analysis(case.system.A)
        dynamic = [m for m in spectrum.modes if abs(m.value) > ZERO_MODE]
        abscissa = float(max(m.real for m in dynamic))
        zeta_min = float(min(m.damping for m in dynamic if m.imag > 0.0))
        reading = read_family(
            anchor, base_case, case, spectrum, common_labels=pinned
        )
        try:
            space = build_action_space(base_case, case, CORE)
            profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
            closure = min(s.margin for s in profile)
        except (np.linalg.LinAlgError, ValueError, InfeasibleReplacement):
            closure = float("nan")
        return {
            "spectral_abscissa": abscissa,
            "zeta_min": zeta_min,
            "alpha_IA": reading.family.alpha if reading else float("nan"),
            "freq_IA_hz": reading.family.frequency_worst_hz if reading else float("nan"),
            "closure_m4": closure,
        }

    def abscissa_of(plan, converter=None):
        case = solve_case(plan, converter=converter)
        spectrum = eigen_analysis(case.system.A)
        return float(
            max(m.real for m in spectrum.modes if abs(m.value) > ZERO_MODE)
        ), case

    unrepaired = measure(flagship)
    baseline = measure(base_case)
    experiment.note(
        f"unrepaired flagship abscissa {unrepaired['spectral_abscissa']:+.4f}, "
        f"closure m4 {unrepaired['closure_m4']:.4f}; "
        f"base closure m4 {baseline['closure_m4']:.4f}"
    )

    rows = []

    def record(strategy, target_name, target, detail, case, **costs):
        payload = measure(case)
        rows.append(
            {
                "strategy": strategy,
                "target": target_name,
                "target_value": target,
                "detail": detail,
                # An optimizer that converges ONTO the constraint lands within its
                # own convergence tolerance of the target, not below it by 1e-9.
                # Judging it at 1e-9 would report a successful repair as a failure.
                "meets_target": bool(payload["spectral_abscissa"] <= target + 1e-6),
                "shortfall": float(payload["spectral_abscissa"] - target),
                "pv_mw_kept": total_pv_mw - costs.get("pv_mw_forgone", 0.0),
                "closure_m4_before": unrepaired["closure_m4"],
                **costs,
                **payload,
            }
        )

    # ---- M1 converter control only -----------------------------------------
    for name, target in zip(TARGET_NAMES, TARGETS, strict=True):
        def constraint(vector, target=target):
            try:
                value, _ = abscissa_of(
                    flagship_plan, converter_from_physical(vector)
                )
                return target - value
            except Exception:  # noqa: BLE001
                return -1.0

        solution = minimize(
            lambda v: float(v @ v),
            np.full(len(PHYSICAL), -0.15),
            method="SLSQP",
            bounds=[(-np.log(8.0), np.log(8.0))] * len(PHYSICAL),
            constraints=[{"type": "ineq", "fun": constraint}],
            options={"maxiter": 90, "ftol": 1e-7},
        )
        converter = converter_from_physical(solution.x)
        case = solve_case(flagship_plan, converter=converter)
        detail = ", ".join(
            f"{n} x{np.exp(v):.3f}"
            for n, v in zip(PHYSICAL, solution.x, strict=True)
            if abs(v) > 0.02
        )
        record(
            "M1 converter only",
            name,
            target,
            detail or "no change",
            case,
            pv_mw_forgone=0.0,
            synchronous_mva_added=0.0,
            controller_change_norm=float(np.linalg.norm(solution.x)),
        )
        print(f"  M1 {name:11s} done", flush=True)

    # ---- M2 condenser only, over location and rating ------------------------
    ratings = np.round(np.arange(0.02, 0.51, 0.02), 4)
    for size in range(1, len(CORE) + 1):
        for where in combinations(CORE, size):
            for rating in ratings:
                try:
                    plan = ReplacementPlan.of(
                        {b: 1.0 for b in CORE}, condenser={b: float(rating) for b in where}
                    )
                    value, case = abscissa_of(plan)
                except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
                    continue
                mva = sum(network.machines[b]["Sn"] for b in where) * float(rating)
                for name, target in zip(TARGET_NAMES, TARGETS, strict=True):
                    if value <= target + 1e-6:
                        record(
                            "M2 condenser only",
                            name,
                            target,
                            f"rating {rating:.2f} at {'+'.join(map(str, where))}",
                            case,
                            pv_mw_forgone=0.0,
                            synchronous_mva_added=float(mva),
                            controller_change_norm=0.0,
                        )
                if value <= min(TARGETS) + 1e-6:
                    break  # smallest rating meeting every target at this location
        print(f"  M2 locations of size {size} done", flush=True)

    # ---- M3 condenser plus minimal retune ----------------------------------
    # A condenser below about a tenth of rating cannot carry the reactive power
    # of the bus at all, so those ratings are infeasible rather than merely
    # insufficient. Feasibility is checked before the retune is optimised on top.
    feasible_ratings = []
    for rating in (0.10, 0.15, 0.20, 0.25):
        try:
            solve_case(
                ReplacementPlan.of(
                    {b: 1.0 for b in CORE}, condenser={b: rating for b in CORE}
                )
            )
            feasible_ratings.append(rating)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            continue
    experiment.note(f"M3 feasible condenser ratings: {feasible_ratings}")
    for rating in feasible_ratings:
        combined_plan = ReplacementPlan.of(
            {b: 1.0 for b in CORE}, condenser={b: rating for b in CORE}
        )
        mva = sum(network.machines[b]["Sn"] for b in CORE) * rating
        for name, target in zip(TARGET_NAMES, TARGETS, strict=True):
            def constraint(vector, target=target, plan=combined_plan):
                try:
                    value, _ = abscissa_of(plan, converter_from_physical(vector))
                    return target - value
                except Exception:  # noqa: BLE001
                    return -1.0

            solution = minimize(
                lambda v: float(v @ v),
                np.full(len(PHYSICAL), -0.10),
                method="SLSQP",
                bounds=[(-np.log(8.0), np.log(8.0))] * len(PHYSICAL),
                constraints=[{"type": "ineq", "fun": constraint}],
                options={"maxiter": 70, "ftol": 1e-7},
            )
            try:
                case = solve_case(
                    combined_plan, converter=converter_from_physical(solution.x)
                )
            except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
                continue
            record(
                "M3 converter + condenser",
                name,
                target,
                f"condenser {rating:.2f} at all four, "
                + (
                    ", ".join(
                        f"{n} x{np.exp(v):.3f}"
                        for n, v in zip(PHYSICAL, solution.x, strict=True)
                        if abs(v) > 0.02
                    )
                    or "no retune"
                ),
                case,
                pv_mw_forgone=0.0,
                synchronous_mva_added=float(mva),
                controller_change_norm=float(np.linalg.norm(solution.x)),
            )
        print(f"  M3 rating {rating} done", flush=True)

    # ---- M4 restore one machine --------------------------------------------
    for bus in CORE:
        plan = ReplacementPlan.of({b: 1.0 for b in CORE if b != bus})
        try:
            value, case = abscissa_of(plan)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            continue
        forgone = float(network.machines[bus]["Sn"] * 1.0)
        for name, target in zip(TARGET_NAMES, TARGETS, strict=True):
            if value <= target + 1e-6:
                record(
                    "M4 restore one machine",
                    name,
                    target,
                    f"machine {bus} kept",
                    case,
                    pv_mw_forgone=forgone,
                    synchronous_mva_added=forgone,
                    controller_change_norm=0.0,
                )
    print("  M4 done", flush=True)

    # ---- M5 fractional synchronous capacity --------------------------------
    for bus in CORE:
        for keep in np.round(np.arange(0.05, 0.96, 0.05), 4):
            mapping = {b: 1.0 for b in CORE}
            mapping[bus] = float(1.0 - keep)
            try:
                value, case = abscissa_of(ReplacementPlan.of(mapping))
            except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
                continue
            forgone = float(network.machines[bus]["Sn"] * keep)
            hit = False
            for name, target in zip(TARGET_NAMES, TARGETS, strict=True):
                if value <= target + 1e-6:
                    record(
                        "M5 fractional machine",
                        name,
                        target,
                        f"keep {keep:.2f} of machine {bus}",
                        case,
                        pv_mw_forgone=forgone,
                        synchronous_mva_added=forgone,
                        controller_change_norm=0.0,
                    )
                    hit = True
            if hit and value <= min(TARGETS) + 1e-6:
                break
    print("  M5 done", flush=True)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E34_mitigation_frontier.csv", parquet=True)

    # ---- Pareto front per target on the three cost axes ---------------------
    def pareto(block):
        axes = block[
            ["pv_mw_forgone", "synchronous_mva_added", "controller_change_norm"]
        ].to_numpy()
        keep = []
        for i in range(len(axes)):
            dominated = np.all(axes <= axes[i], axis=1) & np.any(axes < axes[i], axis=1)
            keep.append(not dominated.any())
        return block[np.array(keep)]

    fronts = []
    for name in TARGET_NAMES:
        block = table[(table.target == name) & table.meets_target]
        if block.empty:
            continue
        front = pareto(block.reset_index(drop=True))
        front = front.assign(target=name)
        fronts.append(front)
    frontier = pd.concat(fronts) if fronts else pd.DataFrame()
    experiment.save_table(frontier, "E34_pareto_front.csv")

    good = table[table.meets_target & table.closure_m4.notna()]
    checks = {
        "candidates": int(len(table)),
        "meeting_target": int(table.meets_target.sum()),
        "pareto_points": int(len(frontier)),
        "closure_m4_unrepaired": unrepaired["closure_m4"],
        "closure_m4_base": baseline["closure_m4"],
        "closure_m4_after_min": float(good.closure_m4.min()),
        "closure_m4_after_median": float(good.closure_m4.median()),
        "repairs_that_increase_closure": int(
            (good.closure_m4 > unrepaired["closure_m4"]).sum()
        ),
        "repairs_total": int(len(good)),
        "strategies_keeping_all_pv": sorted(
            set(table[table.meets_target & (table.pv_mw_forgone == 0)].strategy)
        ),
        "strategies_zero_sync": sorted(
            set(
                table[table.meets_target & (table.synchronous_mva_added == 0)].strategy
            )
        ),
        "total_pv_mw": total_pv_mw,
    }
    verdict = (
        "PASS"
        if checks["repairs_that_increase_closure"] == checks["repairs_total"]
        else "MIXED"
    )
    print()
    for key, value in checks.items():
        print(f"  {key:36s} {value}")
    experiment.finish(verdict, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
