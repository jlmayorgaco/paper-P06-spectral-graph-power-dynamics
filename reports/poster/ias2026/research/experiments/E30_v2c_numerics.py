"""E30 - v2C item 8. Is a small port margin a property of the operator or of the solve?

cond(T_0(j omega)) sits near 9e+03 across the campaign. That is not catastrophic,
but "m_4 becomes small near closure" is the central H5C finding and it must be
shown not to be a linear-solve artifact.

Five checks, declared in the protocol before the run:

  1. backward residual of every T_0 solve
  2. the margin under a doubled frequency grid
  3. the margin under golden-section local minimisation around omega_port
  4. a perturbation ensemble: relative complex perturbation of T_0 at 1e-13,
     giving the spread of m_4 induced by noise at the level of roundoff
  5. three independent solve algorithms, LU, QR and SVD, compared, through an
     INDEPENDENT reimplementation of the closure algebra that is also checked
     against the production one

Extended precision was declared unavailable in advance: mpmath is not installed
in the pinned environment and the interpreter has no pip. Checks 4 and 5 are the
substitute, and the limitation is reported rather than hidden.

The cohorts are those the protocol names: 10 generic accepted samples, the 10
accepted samples nearest the boundary, and the deterministic nominal boundary.

Usage
    python experiments/E30_v2c_numerics.py
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy.linalg import qr, solve_triangular

from _bootstrap import RESULTS
from _v2c_common import BAND_HZ, CORE, base_anchor, nominal_reference, read_family
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import (
    build_action_space,
    closure_at,
    closure_profile,
    refine_closure,
)
from ibr_cycles.uncertainty.sampling import sample_operating_point, stratified_grid

EXPERIMENT = "E30_v2c_numerics"
SEED = 20260912
LOAD_STRATA = np.array([0.90, 0.95, 0.98, 1.00, 1.02, 1.05])
AVAILABILITY_STRATA = np.array([0.75, 0.85, 0.92, 0.96, 1.00])
PER_CELL = 20
BAND_SAMPLES = 61
PERTURBATION = 1e-13
ENSEMBLE = 40


def margin_from_solution(space, s, solution):
    """Closure margin from an already-computed T_0^-1 U, reimplemented here.

    Deliberately independent of ``PortActionSpace.split`` so that agreement
    between the two is evidence about the implementation, not a tautology.
    """

    u = space.selector()
    k = u.T @ solution
    m = space.update(s) @ k
    size = m.shape[0]
    identity = np.eye(size, dtype=np.complex128)
    total = identity + m
    blocks = np.zeros_like(total)
    for a in range(space.order):
        rows = slice(2 * a, 2 * a + 2)
        blocks[rows, rows] = total[rows, rows]
    q = np.linalg.solve(blocks, total) - identity
    return float(np.min(np.abs(np.linalg.eigvals(q) + 1.0)))


def solve_lu(operator, rhs):
    return np.linalg.solve(operator, rhs)


def solve_qr(operator, rhs):
    q, r = qr(operator, mode="economic")
    return solve_triangular(r, q.conj().T @ rhs)


def solve_svd(operator, rhs):
    u, singular, vh = np.linalg.svd(operator, full_matrices=False)
    return vh.conj().T @ ((u.conj().T @ rhs) / singular[:, None])


def audit_sample(space, omega_port, rng, label, extra):
    s = complex(0.0, omega_port)
    operator = space.t0.evaluate(s)
    selector = space.selector()
    production = closure_at(space, omega_port)

    solutions = {
        "lu": solve_lu(operator, selector),
        "qr": solve_qr(operator, selector),
        "svd": solve_svd(operator, selector),
    }
    margins = {
        name: margin_from_solution(space, s, solution)
        for name, solution in solutions.items()
    }
    residuals = {
        name: float(np.linalg.norm(operator @ solution - selector, 2))
        / (float(np.linalg.norm(operator, 2)) * float(np.linalg.norm(solution, 2)))
        for name, solution in solutions.items()
    }

    # one step of iterative refinement on the LU solve
    correction = np.linalg.solve(operator, selector - operator @ solutions["lu"])
    refined_margin = margin_from_solution(space, s, solutions["lu"] + correction)

    # perturbation ensemble at the level of roundoff
    scale = float(np.linalg.norm(operator, "fro")) / np.sqrt(operator.size)
    ensemble = []
    for _ in range(ENSEMBLE):
        noise = PERTURBATION * scale * (
            rng.standard_normal(operator.shape) + 1j * rng.standard_normal(operator.shape)
        )
        try:
            ensemble.append(
                margin_from_solution(space, s, np.linalg.solve(operator + noise, selector))
            )
        except np.linalg.LinAlgError:
            continue
    ensemble = np.array(ensemble)

    coarse = min(closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES),
                 key=lambda x: x.margin)
    fine = min(closure_profile(space, band_hz=BAND_HZ, samples=2 * BAND_SAMPLES - 1),
               key=lambda x: x.margin)
    half_width = 2.0 * np.pi * (BAND_HZ[1] - BAND_HZ[0]) / (BAND_SAMPLES - 1)
    local = refine_closure(space, coarse.omega, half_width=half_width)

    row = {
        "cohort": label,
        "m4_production": production.margin,
        "cond_T0": production.condition,
        "backward_residual": production.solve_residual,
        "m4_lu": margins["lu"],
        "m4_qr": margins["qr"],
        "m4_svd": margins["svd"],
        "algorithm_spread": float(max(margins.values()) - min(margins.values())),
        "reimplementation_gap": abs(margins["lu"] - production.margin),
        "residual_lu": residuals["lu"],
        "residual_qr": residuals["qr"],
        "residual_svd": residuals["svd"],
        "iterative_refinement_shift": abs(refined_margin - margins["lu"]),
        "ensemble_std": float(ensemble.std()) if ensemble.size else float("nan"),
        "ensemble_max_shift": float(np.abs(ensemble - margins["lu"]).max())
        if ensemble.size
        else float("nan"),
        "m4_grid_61": coarse.margin,
        "m4_grid_121": fine.margin,
        "grid_shift": abs(fine.margin - coarse.margin),
        "m4_local_min": local.margin,
        "local_shift": abs(local.margin - coarse.margin),
        "freq_port_grid_hz": coarse.frequency_hz,
        "freq_port_local_hz": local.frequency_hz,
    }
    row.update(extra)
    return row


def main() -> int:
    started = time.time()
    samples = pd.read_csv(
        RESULTS / "tables" / "E28_v2c_modal_family_samples.csv"
    )
    accepted = samples[samples.status == "ACCEPTED"].copy()
    accepted["abs_alpha"] = accepted.alpha_q.abs()
    boundary_ids = set(accepted.nsmallest(10, "abs_alpha")["sample"].astype(int))
    generic_ids = set(
        accepted.sort_values("sample")["sample"].astype(int).to_numpy()[::30][:10]
    )
    wanted = boundary_ids | generic_ids

    base_network = load_network()
    rng = np.random.default_rng(SEED)
    points = stratified_grid(LOAD_STRATA, AVAILABILITY_STRATA, PER_CELL, rng)
    nominal = nominal_reference()
    plan_q = ReplacementPlan.of({b: 1.0 for b in CORE})
    audit_rng = np.random.default_rng(7)

    rows = []
    for index, (load, reactive, availability) in enumerate(points):
        # the sampler must be advanced for EVERY index so the stream matches
        point = sample_operating_point(
            base_network,
            active_load=load,
            reactive_load=reactive,
            availability=availability,
            availability_buses=CORE,
            rng=rng,
        )
        if index not in wanted or not point.dispatchable:
            continue
        try:
            base = solve_case(ReplacementPlan.of({}), network=point.network)
            replaced = solve_case(plan_q, network=point.network)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            continue
        stored = accepted[accepted["sample"] == index]
        if stored.empty:
            continue
        stored = stored.iloc[0]
        space = build_action_space(base, replaced, CORE)
        label = "near boundary" if index in boundary_ids else "generic"
        rows.append(
            audit_sample(
                space,
                float(stored.omega_port),
                audit_rng,
                label,
                {
                    "sample": index,
                    "load": float(stored.load),
                    "availability": float(stored.availability),
                    "alpha_IA": float(stored.alpha_q),
                    "m4_campaign": float(stored.m_4),
                },
            )
        )

    # the deterministic nominal boundary, availability 0.92
    deterministic = np.random.default_rng(0)
    point = sample_operating_point(
        base_network,
        active_load=0.9800,
        reactive_load=1.0,
        availability=0.92,
        availability_buses=CORE,
        rng=deterministic,
        dispatch_jitter=0.0,
        load_scatter=0.0,
    )
    base = solve_case(ReplacementPlan.of({}), network=point.network)
    replaced = solve_case(plan_q, network=point.network)
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)
    reading = read_family(anchor, base, replaced, eigen_analysis(replaced.system.A))
    space = build_action_space(base, replaced, CORE)
    coarse = min(
        closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES),
        key=lambda x: x.margin,
    )
    rows.append(
        audit_sample(
            space,
            coarse.omega,
            audit_rng,
            "deterministic boundary",
            {
                "sample": -1,
                "load": 0.9800,
                "availability": 0.92,
                "alpha_IA": reading.alpha,
                "m4_campaign": float("nan"),
            },
        )
    )

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_audit.csv", index=False)

    summary = {
        "n_audited": len(table),
        "max_backward_residual": float(table.backward_residual.max()),
        "max_algorithm_spread": float(table.algorithm_spread.max()),
        "max_reimplementation_gap": float(table.reimplementation_gap.max()),
        "max_iterative_refinement_shift": float(
            table.iterative_refinement_shift.max()
        ),
        "max_ensemble_std": float(table.ensemble_std.max()),
        "max_ensemble_shift": float(table.ensemble_max_shift.max()),
        "max_grid_shift": float(table.grid_shift.max()),
        "max_local_shift": float(table.local_shift.max()),
        "min_m4_audited": float(table.m4_production.min()),
        "max_cond_T0": float(table.cond_T0.max()),
    }
    headroom = summary["min_m4_audited"] / max(
        summary["max_ensemble_shift"], summary["max_algorithm_spread"], 1e-300
    )

    print("E30 numerical validation of the port margin, %d cases" % len(table))
    print("  cohorts: %s" % dict(table.cohort.value_counts()))
    print()
    print("  smallest m_4 audited            %.6e" % summary["min_m4_audited"])
    print("  worst cond(T0(j omega))         %.3e" % summary["max_cond_T0"])
    print("  worst backward residual         %.3e" % summary["max_backward_residual"])
    print("  LU / QR / SVD max spread        %.3e" % summary["max_algorithm_spread"])
    print("  independent reimplementation    %.3e" % summary["max_reimplementation_gap"])
    print("  one step of refinement moves    %.3e" % summary["max_iterative_refinement_shift"])
    print("  perturbation 1e-13, worst shift %.3e" % summary["max_ensemble_shift"])
    print("  doubled grid moves the margin   %.3e" % summary["max_grid_shift"])
    print("  local minimisation moves it     %.3e" % summary["max_local_shift"])
    print()
    print("  smallest audited margin exceeds the worst numerical shift by a factor"
          " of %.3e" % headroom)
    print()
    print(table[[
        "cohort", "alpha_IA", "m4_production", "cond_T0", "backward_residual",
        "algorithm_spread", "ensemble_max_shift", "grid_shift", "local_shift",
    ]].to_string(index=False, float_format=lambda v: f"{v:11.3e}"))

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=SEED,
        config={
            "perturbation": PERTURBATION,
            "ensemble": ENSEMBLE,
            "band_samples": BAND_SAMPLES,
            "extended_precision": "unavailable: no mpmath and no pip in the pinned venv",
        },
    )
    manifest.finish("REPORTED", headroom_factor=float(headroom), elapsed_s=time.time() - started, **summary)
    print(f"manifest -> {manifest.write(RESULTS / 'manifests')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
