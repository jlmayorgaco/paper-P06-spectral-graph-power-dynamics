"""E21 - PHASE E3. Minimal repair of the Track-A flagship.

Question
    What is the smallest intervention that restores the margin, and how does a
    converter retune compare quantitatively against restoring a synchronous
    machine?

Strategies
    RA  restore one synchronous generator, i.e. undo one replacement
    RB  keep all four replacements, weaken the active-power outer loop only
    RC  keep all four replacements, minimum-norm retune over the converter
        parameters
    RD  keep all four replacements, add the smallest synchronous condenser

Target
    Declared before any optimization: spectral abscissa at most ``-0.05`` over
    non-reference modes. The base case sits at ``-0.126``; the unrepaired
    flagship at ``+0.145``.

Reported for every strategy
    PV megawatts retained, synchronous capacity required, relative parameter
    change, spectral abscissa, tracked inter-area branch, damping, and the
    gauge-invariant closure distance of E17: how far the nearest eigenvalue of Q
    sits from -1.

Status of the result
    NUMERICAL OBSERVATION at one operating point, rho = 1.

Usage
    python experiments/E21_minimal_repair.py
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_admittance import build_action_space

EXPERIMENT = "E21_minimal_repair"
CORE = (30, 33, 35, 37)
TARGET = -0.05
BASE_CONVERTER = ConverterParameters()
TUNED = ("kp_p", "ki_p", "kp_q", "ki_q", "tau_p", "kp_i", "ki_i")


def machine_labels(labels):
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E3 minimal repair")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    reference_spectrum = eigen_analysis(flagship.system.A)
    reference = max(
        (m for m in reference_spectrum.modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    )
    ref_labels = machine_labels(flagship.system.labels)
    ref_index = {n: i for i, n in enumerate(flagship.system.labels)}
    created = reference.value
    total_pv_mw = flagship.replaced_mw

    def evaluate(plan: ReplacementPlan, converter=None) -> dict[str, float]:
        case = solve_case(plan, converter=converter)
        spectrum = eigen_analysis(case.system.A)
        verdict = assess(spectrum, case.system.labels)
        common = sorted(ref_labels & machine_labels(case.system.labels))
        left = reference.right[[ref_index[n] for n in common]]
        index = {n: i for i, n in enumerate(case.system.labels)}
        select = [index[n] for n in common]
        best = max(
            (m for m in spectrum.modes if abs(m.value) > 1e-3),
            key=lambda m: modal_assurance(left, m.right[select]),
        )
        space = build_action_space(base, case, plan.members or CORE)
        try:
            closure = float(abs(space.split(created)["closest_to_minus_one"] + 1.0))
        except np.linalg.LinAlgError:
            closure = float("nan")
        return {
            "spectral_abscissa": verdict.spectral_abscissa,
            "zeta_min": verdict.zeta_min,
            "tracked_real": best.real,
            "tracked_frequency_hz": best.frequency_hz,
            "tracked_mac": modal_assurance(left, best.right[select]),
            "pv_mw_retained": case.replaced_mw,
            "closure_distance": closure,
        }

    def abscissa_of(converter) -> float:
        case = solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE}), converter=converter
        )
        return assess(
            eigen_analysis(case.system.A), case.system.labels
        ).spectral_abscissa

    rows = []

    def record(strategy, detail, payload, change=0.0, sync_mva=0.0):
        rows.append(
            {
                "strategy": strategy,
                "detail": detail,
                "relative_parameter_change": change,
                "synchronous_mva_required": sync_mva,
                "pv_mw_forgone": total_pv_mw - payload["pv_mw_retained"],
                **payload,
            }
        )

    record(
        "R0 none",
        "unrepaired flagship",
        evaluate(ReplacementPlan.of({b: 1.0 for b in CORE})),
    )
    record("R0 none", "base case, all machines", evaluate(ReplacementPlan.of({})))

    # RA restore one synchronous generator
    for bus in CORE:
        members = tuple(b for b in CORE if b != bus)
        payload = evaluate(ReplacementPlan.of({b: 1.0 for b in members}))
        record(
            "RA restore one SG",
            f"machine {bus} kept",
            payload,
            sync_mva=base.dae.network.machines[bus]["Sn"],
        )

    # RB weaken only the active-power outer loop, minimal factor by bisection
    def p_loop(scale: float) -> float:
        return abscissa_of(
            ConverterParameters(
                kp_p=BASE_CONVERTER.kp_p * scale, ki_p=BASE_CONVERTER.ki_p * scale
            )
        )

    scale_star = float(brentq(lambda s: p_loop(s) - TARGET, 0.05, 1.0, xtol=1e-6))
    converter_rb = ConverterParameters(
        kp_p=BASE_CONVERTER.kp_p * scale_star, ki_p=BASE_CONVERTER.ki_p * scale_star
    )
    record(
        "RB weaken P loop",
        f"kp_P and ki_P scaled by {scale_star:.4f}",
        evaluate(ReplacementPlan.of({b: 1.0 for b in CORE}), converter_rb),
        change=float(np.linalg.norm([np.log(scale_star)] * 2)),
    )

    # RC minimum-norm retune over the converter parameters
    baseline = np.array([getattr(BASE_CONVERTER, n) for n in TUNED])

    def build(vector: np.ndarray) -> ConverterParameters:
        values = dict(zip(TUNED, baseline * np.exp(vector), strict=True))
        return ConverterParameters(**values)

    def objective(vector: np.ndarray) -> float:
        return float(vector @ vector)

    def constraint(vector: np.ndarray) -> float:
        try:
            return TARGET - abscissa_of(build(vector))
        except Exception:
            return -1.0

    solution = minimize(
        objective,
        np.full(len(TUNED), -0.2),
        method="SLSQP",
        bounds=[(-np.log(20.0), np.log(20.0))] * len(TUNED),
        constraints=[{"type": "ineq", "fun": constraint}],
        options={"maxiter": 120, "ftol": 1e-8},
    )
    converter_rc = build(solution.x)
    achieved = abscissa_of(converter_rc)
    detail = ", ".join(
        f"{n} x{np.exp(v):.3f}"
        for n, v in zip(TUNED, solution.x, strict=True)
        if abs(v) > 0.02
    )
    record(
        "RC minimum-norm retune",
        detail if achieved <= TARGET + 1e-9 else "target not reached",
        evaluate(ReplacementPlan.of({b: 1.0 for b in CORE}), converter_rc),
        change=float(np.linalg.norm(solution.x)),
    )

    # RD smallest synchronous condenser
    def condenser(rating: float) -> float:
        case = solve_case(
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: rating for b in CORE}
            )
        )
        return assess(
            eigen_analysis(case.system.A), case.system.labels
        ).spectral_abscissa

    ratings = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.25])
    feasible = []
    for rating in ratings:
        try:
            feasible.append((rating, condenser(float(rating))))
        except Exception:
            feasible.append((rating, float("nan")))
    working = [r for r, a in feasible if np.isfinite(a) and a <= TARGET]
    rating_star = float(min(working)) if working else float("nan")
    if np.isfinite(rating_star):
        plan = ReplacementPlan.of(
            {b: 1.0 for b in CORE}, condenser={b: rating_star for b in CORE}
        )
        mva = sum(rating_star * base.dae.network.machines[b]["Sn"] for b in CORE)
        record(
            "RD synchronous condenser",
            f"rating fraction {rating_star:g} at each retired bus",
            evaluate(plan),
            sync_mva=mva,
        )

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_strategies.csv", index=False)
    pd.DataFrame(feasible, columns=["rating", "spectral_abscissa"]).to_csv(
        tables / f"{EXPERIMENT}_condenser_sweep.csv", index=False
    )

    # The bisected strategies land exactly on the target, so the acceptance test
    # carries a tolerance; without it a root-finder success reads as a failure.
    succeeded = table[
        (table.spectral_abscissa <= TARGET + 1e-6) & (table.strategy != "R0 none")
    ]
    keeps_all = succeeded[succeeded.pv_mw_forgone < 1.0]

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"core": list(CORE), "target_abscissa": TARGET, "tuned": list(TUNED)},
    )
    manifest.finish(
        "SUCCESS",
        total_pv_mw=total_pv_mw,
        minimal_p_loop_scale=scale_star,
        minimum_norm_change=float(np.linalg.norm(solution.x)),
        minimal_condenser_rating=rating_star,
        strategies=table.to_dict("records"),
        strategies_meeting_target=len(succeeded),
        strategies_keeping_all_pv=len(keeps_all),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    columns = [
        "strategy",
        "detail",
        "pv_mw_retained",
        "pv_mw_forgone",
        "synchronous_mva_required",
        "relative_parameter_change",
        "spectral_abscissa",
        "tracked_real",
        "zeta_min",
        "closure_distance",
    ]
    print(
        "target: spectral abscissa <= %+.3f   total PV replaced: %.1f MW"
        % (TARGET, total_pv_mw)
    )
    print()
    print(table[columns].to_string(index=False, float_format=lambda v: f"{v:11.4f}"))
    print()
    print("  minimal P-loop weakening that reaches the target : x%.4f" % scale_star)
    print(
        "  minimum-norm retune, relative log change         : %.4f"
        % np.linalg.norm(solution.x)
    )
    print(
        "  smallest condenser rating that reaches the target: %s"
        % (f"{rating_star:g}" if np.isfinite(rating_star) else "none in the sweep")
    )
    print("  strategies meeting the target                    : %d" % len(succeeded))
    print("  of those, keeping every PV megawatt              : %d" % len(keeps_all))
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
