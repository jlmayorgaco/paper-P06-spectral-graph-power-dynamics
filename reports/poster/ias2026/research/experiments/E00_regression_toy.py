"""E00 - GATE 0.

Question
    Does the framework reproduce the frozen five-state regression envelope, and
    do the algebraic identities it relies on hold to machine precision?

Status of the result
    EXACT IDENTITY for the Schur, determinant-lemma and factorization checks.
    NUMERICAL OBSERVATION for the spectral envelope, the Moebius orders and the
    winding provenance, all of which are properties of the frozen amplitudes.

Usage
    python experiments/E00_regression_toy.py
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.actions.interaction import split_interaction
from ibr_cycles.actions.portfolios import (
    PortfolioOutcome,
    minimum_destabilizing_order,
    subsets,
)
from ibr_cycles.cycles.cumulants import convergence_study
from ibr_cycles.cycles.graph import build_action_graph, rank_cycles
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.cycles.winding import (
    check_admissibility,
    circle,
    right_half_plane,
    spectra_union,
    winding_number,
)
from ibr_cycles.diagnosis.provenance import decompose_winding
from ibr_cycles.diagnosis.repair import scalar_repair
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models import toy5
from ibr_cycles.models.toy5_case import (
    ABSCISSA_TOLERANCE,
    ACTION_NAMES,
    EXPECTED_ABSCISSA,
    EXPECTED_K_CRITICAL,
    EXPECTED_K_REPAIR,
    PROBE_POINTS,
    build_case,
)
from ibr_cycles.reduction.residues import expand_self_energy
from ibr_cycles.reduction.schur import schur_determinant_residual
from ibr_cycles.reduction.self_energy import SelfEnergy

EXPERIMENT = "E00_regression_toy"


def portfolio_table(case) -> pd.DataFrame:
    rows = []
    for members in subsets(ACTION_NAMES):
        spectrum = eigen_analysis(case.state_matrix(members))
        critical = spectrum.critical
        expected = EXPECTED_ABSCISSA[members]
        rows.append(
            {
                "portfolio": "".join(members) or "BASE",
                "size": len(members),
                "spectral_abscissa": spectrum.spectral_abscissa,
                "expected_abscissa": expected,
                "deviation": spectrum.spectral_abscissa - expected,
                "critical_imag": critical.imag,
                "frequency_hz": critical.frequency_hz,
                "damping": critical.damping,
                "rhp_modes": spectrum.rhp_count(),
                "eigenvalue_condition": critical.condition,
                "stable": spectrum.is_stable,
            }
        )
    return pd.DataFrame(rows)


def outcomes(case) -> list[PortfolioOutcome]:
    result = []
    for members in subsets(ACTION_NAMES, include_empty=False):
        spectrum = eigen_analysis(case.state_matrix(members))
        critical = spectrum.critical
        result.append(
            PortfolioOutcome(
                members=members,
                spectral_abscissa=spectrum.spectral_abscissa,
                critical_real=critical.real,
                critical_imag=critical.imag,
                damping=critical.damping,
                frequency_hz=critical.frequency_hz,
                rhp_modes=spectrum.rhp_count(),
                eigenvalue_condition=critical.condition,
            )
        )
    return result


def identity_checks(case) -> dict[str, float]:
    green = case.green(ACTION_NAMES)
    sigma = SelfEnergy.from_matrix(case.a0, case.partition)
    expansion = expand_self_energy(sigma)
    schur = max(
        schur_determinant_residual(
            case.state_matrix(m), case.partition, PROBE_POINTS
        ).max_relative_error
        for m in subsets(ACTION_NAMES)
    )
    angle = 0.0
    for s in PROBE_POINTS:
        full = np.linalg.det(s * np.eye(5) - case.a0)
        factored = (s + case.parameters.a) * np.linalg.det(
            toy5.angle_operator(case.parameters, s)
        )
        angle = max(angle, abs(full - factored) / max(abs(full), 1e-300))
    factorization = max(
        split_interaction(green, s).factorization_error for s in PROBE_POINTS
    )
    return {
        "schur_determinant_max_relative_error": schur,
        "angle_self_energy_max_relative_error": angle,
        "determinant_lemma_max_absolute_error": green.lemma_residual(PROBE_POINTS),
        "individual_collective_factorization_error": factorization,
        "self_energy_pole_expansion_error": expansion.max_reconstruction_error,
        "max_action_reconstruction_error": max(
            f.reconstruction_error for f in case.factors
        ),
        "max_action_numerical_rank": max(f.numerical_rank for f in case.factors),
    }


def provenance_study(case, created: complex) -> dict[str, object]:
    green = case.green(ACTION_NAMES)
    lower = spectra_union(
        [
            eigen_analysis(case.state_matrix(m)).values
            for m in subsets(ACTION_NAMES, include_empty=False)
            if len(m) < 3
        ]
    )
    base = eigen_analysis(case.state_matrix(())).values
    portfolio = eigen_analysis(case.state_matrix(ACTION_NAMES)).values

    rows = []
    for radius in (0.10, 0.15, 0.20, 0.30):
        contour = circle(created, radius, samples=2048)
        admissibility = check_admissibility(
            contour,
            portfolio_modes=portfolio,
            base_modes=base,
            lower_order_modes=lower,
        )
        provenance = decompose_winding(green, contour, admissibility)
        rows.append({"region": f"circle_r{radius:g}", **provenance.as_row()})

    half_plane = right_half_plane(boundary=-5e-3, radius=60.0, samples=8192)
    subset_windings = {}
    for members in subsets(ACTION_NAMES, include_empty=False):
        sub = case.green(members)
        result = winding_number(lambda s, g=sub: g.determinant_ratio(s), half_plane)
        subset_windings["".join(members)] = result.nearest_integer
    return {"regions": rows, "right_half_plane_windings": subset_windings}


def main() -> int:
    parser = argparse.ArgumentParser(description="GATE 0 five-state regression")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()

    case = build_case()
    config = {
        "model": "toy5",
        "parameters": dict(case.parameters.__dict__),
        "actions": toy5.ACTIONS,
        "probe_points": [str(s) for s in PROBE_POINTS],
        "abscissa_tolerance": ABSCISSA_TOLERANCE,
    }
    manifest = Manifest(experiment=EXPERIMENT, seed=args.seed, config=config)

    table = portfolio_table(case)
    identities = identity_checks(case)
    green = case.green(ACTION_NAMES)
    row = table.loc[table.portfolio == "ABC"].iloc[0]
    created = complex(float(row.spectral_abscissa), float(row.critical_imag))

    graph = build_action_graph(green, created)
    cycles = rank_cycles(green, created, graph.simple_cycles())
    mobius = decompose(
        lambda m: eigen_analysis(case.state_matrix(tuple(m))).spectral_abscissa,
        ACTION_NAMES,
    )
    cumulants = {
        "pair_AB": convergence_study(green, (0, 1), created).best_error,
        "pair_AC": convergence_study(green, (0, 2), created).best_error,
        "pair_BC": convergence_study(green, (1, 2), created).best_error,
        "triple_ABC": convergence_study(green, (0, 1, 2), created).best_error,
    }

    stressed = toy5.apply_actions(case.parameters, ("A", "B"))
    critical_k, repaired_k = scalar_repair(
        lambda k: (
            eigen_analysis(toy5.state_matrix(replace(stressed, k=k))).spectral_abscissa
        ),
        (0.8, 1.57),
        target=-0.02,
    )

    provenance = provenance_study(case, created)
    max_deviation = float(table.deviation.abs().max())
    kappa = minimum_destabilizing_order(outcomes(case))

    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_portfolios.csv", index=False)
    pd.DataFrame(provenance["regions"]).to_csv(
        tables / f"{EXPERIMENT}_provenance.csv", index=False
    )
    ordered = sorted(mobius.terms.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))
    pd.DataFrame(
        [
            {
                "order": len(key),
                "members": "".join(sorted(key)) or "BASE",
                "term": value,
            }
            for key, value in ordered
        ]
    ).to_csv(tables / f"{EXPERIMENT}_mobius.csv", index=False)

    gate_passed = bool(
        max_deviation < ABSCISSA_TOLERANCE
        and kappa == 3
        and identities["schur_determinant_max_relative_error"] < 1e-10
        and identities["determinant_lemma_max_absolute_error"] < 1e-8
        and abs(critical_k - EXPECTED_K_CRITICAL) < 1e-6
        and abs(repaired_k - EXPECTED_K_REPAIR) < 1e-6
    )

    manifest.finish(
        "SUCCESS" if gate_passed else "GATE_FAILED",
        max_abscissa_deviation=max_deviation,
        minimum_destabilizing_order=kappa,
        identities=identities,
        cumulant_best_relative_errors=cumulants,
        mobius_third_order_term=mobius.terms[frozenset(ACTION_NAMES)],
        mobius_truncation_order_1=mobius.truncated(ACTION_NAMES, 1),
        mobius_truncation_order_2=mobius.truncated(ACTION_NAMES, 2),
        mobius_truncation_order_3=mobius.truncated(ACTION_NAMES, 3),
        action_graph_feedback_cores=graph.feedback_cores,
        strongest_cycle=(
            {"members": cycles[0].members, "magnitude": cycles[0].magnitude}
            if cycles
            else None
        ),
        controller_surgery={
            "k_critical": critical_k,
            "k_repair_target_minus_0_02": repaired_k,
            "k_unsafe": case.parameters.k + toy5.ACTIONS["C"]["k"],
        },
        provenance=provenance,
    )
    path = manifest.write(RESULTS / "manifests")

    print(table.to_string(index=False))
    print()
    print(json.dumps(identities, indent=2))
    print()
    print(f"kappa = {kappa}   max abscissa deviation = {max_deviation:.3e}")
    print(
        "Moebius truncations: order1 = "
        f"{mobius.truncated(ACTION_NAMES, 1):+.6f}, order2 = "
        f"{mobius.truncated(ACTION_NAMES, 2):+.6f}, exact = "
        f"{mobius.truncated(ACTION_NAMES, 3):+.6f}"
    )
    print(f"controller surgery: k_crit = {critical_k:.6f}, k_repair = {repaired_k:.6f}")
    print(f"manifest -> {path}")
    print("GATE 0 PASSED" if gate_passed else "GATE 0 FAILED")
    return 0 if gate_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
