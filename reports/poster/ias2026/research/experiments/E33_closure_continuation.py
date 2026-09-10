"""E33 - overnight 5. Dense closure continuation against the full-system boundary.

Deterministic continuation over load and photovoltaic availability, using the
frozen v2C modal-family observable and the frozen 61-point closure metric. The
question is whether the port-space closure sees the same boundary the full
eigenvalue problem sees:

    alpha_IA -> 0   while   m4 -> 0        and    omega_port/(2 pi) ~= f_IA

The frozen metric is used for the verdict. Local refinement of the band minimum
is computed alongside as a SECONDARY numerical validation, because E30 showed the
61-point grid is a conservative upper bound near the boundary. The frozen metric
is not replaced.
"""

from __future__ import annotations

import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, choose_workers, pin_blas_threads
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
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import (
    build_action_space,
    closure_profile,
    refine_closure,
)
from ibr_cycles.uncertainty.sampling import sample_operating_point

NAME = "E33_closure_continuation"
BAND_SAMPLES = 61
LOADS = np.round(np.arange(0.900, 1.0501, 0.0025), 6)
AVAILABILITIES = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)

_CONTEXT: dict = {}


def _initialise():
    pin_blas_threads()
    base_network = load_network()
    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    _CONTEXT.update(
        network=base_network,
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
        plan=ReplacementPlan.of({b: 1.0 for b in CORE}),
    )


def _evaluate(task):
    load, availability = task
    network = _CONTEXT["network"]
    row = {"load": float(load), "availability": float(availability)}
    rng = np.random.default_rng(0)  # deterministic: no jitter, no scatter
    point = sample_operating_point(
        network,
        active_load=float(load),
        reactive_load=1.0,
        availability=float(availability),
        availability_buses=CORE,
        rng=rng,
        dispatch_jitter=0.0,
        load_scatter=0.0,
    )
    row["dispatchable"] = bool(point.dispatchable)
    if not point.dispatchable:
        row["status"] = "NOT_DISPATCHABLE"
        return row
    try:
        base = solve_case(ReplacementPlan.of({}), network=point.network)
        replaced = solve_case(_CONTEXT["plan"], network=point.network)
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
        row["status"] = f"INFEASIBLE: {str(error)[:70]}"
        return row
    spectrum = eigen_analysis(replaced.system.A)
    anchor = base_anchor(eigen_analysis(base.system.A), base, _CONTEXT["nominal"])
    if anchor is None:
        row["status"] = "TRACKING_FAILURE"
        return row
    reading = read_family(
        anchor, base, replaced, spectrum, common_labels=_CONTEXT["pinned"]
    )
    if reading is None:
        row["status"] = "TRACKING_FAILURE"
        return row
    family = reading.family
    row.update(
        {
            "status": "OK",
            "alpha_IA": family.alpha,
            "freq_IA_hz": family.frequency_worst_hz,
            "alpha_base": float(anchor.real),
            "family_size": family.size,
        }
    )
    try:
        space = build_action_space(base, replaced, CORE)
        profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
        best = min(profile, key=lambda s: s.margin)
        half_width = 2.0 * np.pi * (BAND_HZ[1] - BAND_HZ[0]) / (BAND_SAMPLES - 1)
        refined = refine_closure(space, best.omega, half_width=half_width)
        split = space.split(complex(0.0, best.omega))
        mu = complex(split["closest_to_minus_one"])
        row.update(
            {
                "m4_frozen": best.margin,
                "omega_port": best.omega,
                "freq_port_hz": best.frequency_hz,
                "cond_T0_at_omega_port": best.condition,
                "solve_residual": best.solve_residual,
                "worst_cond_T0": max(s.condition for s in profile),
                "m4_refined": refined.margin,
                "freq_port_refined_hz": refined.frequency_hz,
                "mu_real": mu.real,
                "mu_imag": mu.imag,
            }
        )
    except (np.linalg.LinAlgError, ValueError) as error:
        row["closure_error"] = str(error)[:70]
    return row


def main() -> int:
    pin_blas_threads()
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Does the port closure see the same boundary the full eigenproblem sees?",
        config={
            "loads": [float(LOADS.min()), float(LOADS.max()), 0.0025],
            "availabilities": list(AVAILABILITIES),
            "band_samples_frozen": BAND_SAMPLES,
            "refinement": "golden section, SECONDARY only, frozen metric not replaced",
            "deterministic": True,
        },
        workers=workers,
    )
    started = time.time()
    tasks = [(load, availability) for availability in AVAILABILITIES for load in LOADS]
    experiment.note(f"{len(tasks)} continuation points on {workers} workers")

    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, tasks, chunksize=4)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E33_closure_continuation.csv", parquet=True)
    good = table[(table.status == "OK") & table.m4_frozen.notna()].copy()
    good["abs_alpha"] = good.alpha_IA.abs()
    good["delta_f_hz"] = good.freq_port_hz - good.freq_IA_hz
    good["delta_f_refined_hz"] = good.freq_port_refined_hz - good.freq_IA_hz

    frozen_rho = spearmanr(good.m4_frozen, good.abs_alpha)
    refined_rho = spearmanr(good.m4_refined, good.abs_alpha)
    near = good[good.abs_alpha <= 0.05]
    slope, intercept = np.polyfit(good.abs_alpha, np.log(good.m4_frozen), 1)

    # boundary error: where does each curve cross zero, per availability
    boundary_rows = []
    for availability in AVAILABILITIES:
        block = good[good.availability == availability].sort_values("load")
        if block.empty:
            continue
        sign = np.sign(block.alpha_IA.to_numpy())
        crossings = np.flatnonzero(np.diff(sign) > 0)
        alpha_zero = (
            float(block.load.to_numpy()[crossings[0]]) if crossings.size else np.nan
        )
        m4_min_load = float(block.load.to_numpy()[int(np.argmin(block.m4_frozen))])
        m4_refined_min_load = float(
            block.load.to_numpy()[int(np.argmin(block.m4_refined))]
        )
        boundary_rows.append(
            {
                "availability": availability,
                "load_alpha_zero": alpha_zero,
                "load_m4_minimum_frozen": m4_min_load,
                "load_m4_minimum_refined": m4_refined_min_load,
                "boundary_error_frozen": abs(m4_min_load - alpha_zero),
                "boundary_error_refined": abs(m4_refined_min_load - alpha_zero),
                "m4_at_boundary_frozen": float(
                    block.iloc[crossings[0]].m4_frozen if crossings.size else np.nan
                ),
                "m4_at_boundary_refined": float(
                    block.iloc[crossings[0]].m4_refined if crossings.size else np.nan
                ),
                "points": int(len(block)),
            }
        )
    boundaries = pd.DataFrame(boundary_rows)
    experiment.save_table(boundaries, "E33_boundaries.csv")

    checks = {
        "points_total": int(len(table)),
        "points_dispatchable": int(table.dispatchable.sum()),
        "points_evaluated": int(len(good)),
        "spearman_m4_frozen_vs_abs_alpha": float(frozen_rho.statistic),
        "spearman_m4_frozen_p": float(frozen_rho.pvalue),
        "spearman_m4_refined_vs_abs_alpha": float(refined_rho.statistic),
        "log_m4_slope_on_abs_alpha": float(slope),
        "median_abs_delta_f_hz": float(good.delta_f_hz.abs().median()),
        "p95_abs_delta_f_hz": float(np.percentile(good.delta_f_hz.abs(), 95)),
        "median_abs_delta_f_near_boundary_hz": float(near.delta_f_hz.abs().median()),
        "p95_abs_delta_f_near_boundary_hz": float(
            np.percentile(near.delta_f_hz.abs(), 95)
        ),
        "n_near_boundary": int(len(near)),
        "min_m4_frozen_near_boundary": float(near.m4_frozen.min()),
        "min_m4_refined_near_boundary": float(near.m4_refined.min()),
        "median_boundary_error_frozen": float(boundaries.boundary_error_frozen.median()),
        "median_boundary_error_refined": float(
            boundaries.boundary_error_refined.median()
        ),
        "worst_cond_T0": float(good.worst_cond_T0.max()),
        "worst_solve_residual": float(good.solve_residual.max()),
        "mu_closest_approach_to_minus_one": float(
            np.min(np.abs((good.mu_real + 1j * good.mu_imag) + 1.0))
        ),
    }
    verdict = (
        "PASS"
        if checks["spearman_m4_frozen_vs_abs_alpha"] > 0.8
        and checks["median_abs_delta_f_near_boundary_hz"] < 0.05
        else "MIXED"
    )
    print()
    for key, value in checks.items():
        print(f"  {key:42s} {value}")
    print()
    print(boundaries.to_string(index=False, float_format=lambda v: f"{v:9.5f}"))
    experiment.finish(verdict, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
