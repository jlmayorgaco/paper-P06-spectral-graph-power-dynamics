"""E01 - GATE 1.

Question
    On a realistic converter-dominated case, does every single action and every
    pair stay stable while the full portfolio does not, and if so which factor
    of the determinant carries the created mode?

Hypothesis under test
    That the created mode belongs to the collective interaction factor rather
    than to any isolated action.

Independent variables
    the eight subsets of the three stress actions A (SG damping to zero),
    B (active-power loop gain 0.20 -> 2.00), C (PLL retuned to wn = 12 rad/s,
    zeta = 0.15), at a single frozen AC operating point.

Dependent variables
    spectral abscissa, critical frequency and damping, winding of the full,
    individual and collective determinant factors, and the Moebius interaction
    orders of the spectral abscissa.

Controls
    equilibrium residuals certified per portfolio; additivity of the reduced
    matrix measured; operating-point drift measured; every identity checked
    against a direct determinant.

Status of the result
    EXACT IDENTITY for the algebraic checks. NUMERICAL OBSERVATION for the
    envelope, the provenance and the interaction orders, at these amplitudes
    and this operating point only.

Usage
    python experiments/E01_regression_ieee9.py [--skip-halfplane]
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.actions.interaction import split_interaction
from ibr_cycles.actions.portfolios import (
    PortfolioOutcome,
    minimum_destabilizing_order,
    subsets,
)
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
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models import ieee9_gfl
from ibr_cycles.models.ieee9_case import (
    ABSCISSA_TOLERANCE,
    ACTION_NAMES,
    EXPECTED_ABSCISSA,
    PROBE_POINTS,
    RHP_BOUNDARY,
    RHP_RADIUS,
    SPECIFICATION_ABSCISSA,
    SPECIFICATION_TOLERANCE,
    build_case,
)
from ibr_cycles.reduction.residues import expand_self_energy
from ibr_cycles.reduction.schur import schur_determinant_residual
from ibr_cycles.reduction.self_energy import SelfEnergy

EXPERIMENT = "E01_regression_ieee9"


def portfolio_table(case) -> pd.DataFrame:
    rows = []
    for members in subsets(ACTION_NAMES):
        equilibrium = case.equilibrium(members)
        spectrum = eigen_analysis(case.state_matrix(members))
        critical = spectrum.critical
        rows.append(
            {
                "portfolio": "".join(members) or "BASE",
                "size": len(members),
                "spectral_abscissa": spectrum.spectral_abscissa,
                "frozen_expected": EXPECTED_ABSCISSA[members],
                "specification_reference": SPECIFICATION_ABSCISSA[members],
                "deviation_frozen": spectrum.spectral_abscissa
                - EXPECTED_ABSCISSA[members],
                "deviation_specification": spectrum.spectral_abscissa
                - SPECIFICATION_ABSCISSA[members],
                "critical_imag": critical.imag,
                "frequency_hz": critical.frequency_hz,
                "damping": critical.damping,
                "rhp_modes": spectrum.rhp_count(),
                "eigenvalue_condition": critical.condition,
                "equilibrium_norm_f": equilibrium.norm_f,
                "equilibrium_norm_g": equilibrium.norm_g,
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
    green = case.green()
    sigma = SelfEnergy.from_matrix(case.a0, case.partition)
    expansion = expand_self_energy(sigma)
    critical = eigen_analysis(case.state_matrix(ACTION_NAMES)).critical.value
    split = split_interaction(green, critical)
    return {
        "self_energy_pole_expansion_status": expansion.status,
        "self_energy_hidden_eigenvalue_condition": expansion.poles[0].condition,
        "schur_determinant_max_relative_error": max(
            schur_determinant_residual(
                case.state_matrix(m), case.partition, PROBE_POINTS
            ).max_relative_error
            for m in subsets(ACTION_NAMES)
        ),
        "determinant_lemma_max_absolute_error": green.lemma_residual(PROBE_POINTS),
        "individual_collective_factorization_error": max(
            split_interaction(green, s).factorization_error for s in PROBE_POINTS
        ),
        "reduced_matrix_additivity_residual": case.additivity_residual(),
        "operating_point_drift": case.equilibrium_drift(),
        "max_action_reconstruction_error": max(
            f.reconstruction_error for f in case.factors()
        ),
        "max_action_numerical_rank": max(f.numerical_rank for f in case.factors()),
        "abs_full_at_created_mode": abs(split.full),
        "abs_individual_at_created_mode": abs(split.individual),
        "abs_collective_at_created_mode": abs(split.collective),
        "min_singular_value_of_I_plus_Q": split.min_singular_value,
        "spectral_radius_of_Q": split.spectral_radius_q,
    }


def provenance_rows(case) -> list[dict[str, object]]:
    green = case.green()
    created = eigen_analysis(case.state_matrix(ACTION_NAMES)).critical.value
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
    for radius in (0.005, 0.010, 0.020, 0.050):
        contour = circle(created, radius, samples=4096)
        admissibility = check_admissibility(
            contour,
            portfolio_modes=portfolio,
            base_modes=base,
            lower_order_modes=lower,
        )
        result = decompose_winding(green, contour, admissibility)
        rows.append({"region": f"circle_r{radius:g}", **result.as_row()})
    return rows


def half_plane_windings(case) -> dict[str, object]:
    """Argument-principle cross-check of the stability verdict on Re(s) > b.

    The sample step is tied to the distance from the contour to the nearest
    mode. Without that tie every portfolio returns a clean integer zero, because
    the phase excursion around a lightly damped mode cancels between samples.
    """

    green = case.green()
    modes = np.concatenate(
        [eigen_analysis(case.state_matrix(m)).values for m in subsets(ACTION_NAMES)]
    )
    clearance = float(np.abs(modes.real - RHP_BOUNDARY).min())
    contour = right_half_plane(boundary=RHP_BOUNDARY, radius=RHP_RADIUS, samples=1200)
    counts: dict[str, int] = {}
    resolved = True
    for members in subsets(ACTION_NAMES, include_empty=False):
        subset = green.subset(members)
        result = winding_number(
            lambda s, g=subset: g.determinant_ratio(s),
            contour,
            max_segment=clearance / 4.0,
            refine=False,
        )
        counts["".join(members)] = result.nearest_integer
        resolved = resolved and result.resolved_for(clearance)
    return {
        "boundary": RHP_BOUNDARY,
        "clearance": clearance,
        "max_segment": clearance / 4.0,
        "resolved": resolved,
        "counts": counts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="GATE 1 IEEE-9 GFL regression")
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument(
        "--skip-halfplane",
        action="store_true",
        help="skip the resolved right-half-plane winding cross-check (slow)",
    )
    args = parser.parse_args()

    case = build_case()
    config = {
        "model": "ieee9_gfl",
        "parameters": dict(case.parameters.__dict__),
        "actions": ieee9_gfl.ACTIONS,
        "loads": {str(k): str(v) for k, v in ieee9_gfl.LOADS.items()},
        "retained_states": ieee9_gfl.RETAINED_STATES,
        "probe_points": [str(s) for s in PROBE_POINTS],
    }
    manifest = Manifest(experiment=EXPERIMENT, seed=args.seed, config=config)

    table = portfolio_table(case)
    identities = identity_checks(case)
    provenance = provenance_rows(case)
    kappa = minimum_destabilizing_order(outcomes(case))
    green = case.green()
    created = eigen_analysis(case.state_matrix(ACTION_NAMES)).critical.value

    graph = build_action_graph(green, created)
    cycles = rank_cycles(green, created, graph.simple_cycles())
    mobius = decompose(
        lambda m: eigen_analysis(case.state_matrix(tuple(m))).spectral_abscissa,
        ACTION_NAMES,
    )
    truncations = {
        "order_1": mobius.truncated(ACTION_NAMES, 1),
        "order_2": mobius.truncated(ACTION_NAMES, 2),
        "exact": mobius.truncated(ACTION_NAMES, 3),
        "irreducible_third_order": mobius.terms[frozenset(ACTION_NAMES)],
    }
    mechanism = (
        "GENUINE_THIRD_ORDER"
        if truncations["order_2"] < 0.0 < truncations["exact"]
        else "PAIRWISE_ACCUMULATION"
        if truncations["order_1"] < 0.0 < truncations["order_2"]
        else "SINGLE_ACTION_DOMINATED"
    )

    half_plane = None if args.skip_halfplane else half_plane_windings(case)

    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_portfolios.csv", index=False)
    pd.DataFrame(provenance).to_csv(
        tables / f"{EXPERIMENT}_provenance.csv", index=False
    )
    ordered = sorted(mobius.terms.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))
    pd.DataFrame(
        [
            {"order": len(k), "members": "".join(sorted(k)) or "BASE", "term": v}
            for k, v in ordered
        ]
    ).to_csv(tables / f"{EXPERIMENT}_mobius.csv", index=False)

    frozen_ok = float(table.deviation_frozen.abs().max()) < ABSCISSA_TOLERANCE
    spec_ok = float(table.deviation_specification.abs().max()) < SPECIFICATION_TOLERANCE
    collective = [row for row in provenance if row["verdict"] == "COLLECTIVE"]
    gate_passed = bool(
        frozen_ok
        and kappa == 3
        and identities["schur_determinant_max_relative_error"] < 1e-10
        and identities["determinant_lemma_max_absolute_error"] < 1e-10
        and identities["reduced_matrix_additivity_residual"] < 1e-12
        and len(collective) >= 2
    )

    manifest.finish(
        "SUCCESS" if gate_passed else "GATE_FAILED",
        max_deviation_from_frozen=float(table.deviation_frozen.abs().max()),
        max_deviation_from_specification=float(
            table.deviation_specification.abs().max()
        ),
        specification_within_tolerance=spec_ok,
        minimum_destabilizing_order=kappa,
        identities=identities,
        provenance=provenance,
        right_half_plane=half_plane,
        mobius_terms={"".join(sorted(k)) or "BASE": v for k, v in ordered},
        mobius_truncations=truncations,
        mechanism_verdict=mechanism,
        action_graph_feedback_cores=graph.feedback_cores,
        strongest_cycles=[
            {"members": c.members, "magnitude": c.magnitude} for c in cycles[:3]
        ],
        power_flow={
            "max_mismatch": case.power_flow.max_mismatch,
            "min_voltage": case.power_flow.min_voltage,
            "max_voltage": case.power_flow.max_voltage,
        },
    )
    path = manifest.write(RESULTS / "manifests")

    columns = [
        "portfolio",
        "spectral_abscissa",
        "deviation_frozen",
        "deviation_specification",
        "critical_imag",
        "frequency_hz",
        "damping",
        "stable",
    ]
    print(table[columns].to_string(index=False))
    print()
    print(json.dumps(identities, indent=2))
    print()
    print(f"kappa = {kappa}")
    print(
        "Moebius truncations: order1 = {order_1:+.8f}, order2 = {order_2:+.8f}, "
        "exact = {exact:+.8f}, irreducible third order = "
        "{irreducible_third_order:+.8f}".format(**truncations)
    )
    print(f"mechanism verdict: {mechanism}")
    for row in provenance:
        print(
            f"  {row['region']:>14}: {row['verdict']:>20}  "
            f"full={row['integer_full']:+d} individual={row['integer_individual']:+d} "
            f"collective={row['integer_collective']:+d}"
        )
    if half_plane is not None:
        boundary = half_plane["boundary"]
        print(f"  right half plane Re(s) > {boundary}: {half_plane['counts']}")
    print(f"manifest -> {path}")
    print("GATE 1 PASSED" if gate_passed else "GATE 1 FAILED")
    return 0 if gate_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
