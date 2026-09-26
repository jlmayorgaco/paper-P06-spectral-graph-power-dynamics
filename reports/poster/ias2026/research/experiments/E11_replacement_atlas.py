"""E11 - PHASE B. Individual SG to PV-GFL replacement atlas (R1).

Question
    How much of each synchronous machine can be replaced by a battery-free
    PV grid-following converter before the dynamics degrade, and does a nodal
    strength ranking predict that answer?

Independent variables
    candidate bus (9), replacement fraction rho in {0.25, 0.50, 0.75, 1.00},
    reactive policy in {matched, unity_pf}.

Dependent variables
    feasibility, device loading, spectral abscissa over non-reference modes,
    minimum damping in 0.1-5 Hz, critical frequency, mode condition number.

Controls
    equilibrium residual and index-1 conditioning per case; operating-point
    drift measured for the matched policy, where it must be exactly zero.

Status of the result
    NUMERICAL OBSERVATION at one operating point and one converter tuning.

Usage
    python experiments/E11_replacement_atlas.py
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.baselines import nodal_metrics
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

EXPERIMENT = "E11_replacement_atlas"
RHO_GRID = (0.25, 0.50, 0.75, 1.00)
POLICIES = ("matched", "unity_pf")
MATERIAL_DEGRADATION = 0.005


def evaluate(plan: ReplacementPlan, base_z0: np.ndarray | None) -> dict[str, object]:
    case = solve_case(plan)
    verdict = assess(eigen_analysis(case.system.A), case.system.labels)
    drift = (
        float(np.abs(case.equilibrium.z - base_z0).max())
        if base_z0 is not None
        else float("nan")
    )
    gfl_loading = max(
        (slot.loading for slot in case.dae.slots if slot.kind == "gfl"), default=0.0
    )
    sg_loading = max(
        (slot.loading for slot in case.dae.slots if slot.kind == "sg"), default=0.0
    )
    return {
        "n_states": case.n_states,
        "gz_condition": case.gz_condition,
        "equilibrium_norm_f": case.equilibrium.norm_f,
        "equilibrium_norm_g": case.equilibrium.norm_g,
        "operating_point_drift": drift,
        "gfl_loading": gfl_loading,
        "sg_loading": sg_loading,
        "replaced_mw": case.replaced_mw,
        "spectral_abscissa": verdict.spectral_abscissa,
        "zeta_min": verdict.zeta_min,
        "zeta_min_frequency_hz": verdict.zeta_min_frequency_hz,
        "critical_frequency_hz": verdict.critical_frequency_hz,
        "critical_condition": verdict.critical_condition,
        "rhp_modes": verdict.rhp_modes,
        "unstable": verdict.unstable,
        "reference_clean": verdict.reference_exclusion_is_clean,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE B replacement atlas")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    network = load_network()
    power_flow = solve_power_flow(network)
    metrics = nodal_metrics(network, power_flow)
    base_case = solve_case(ReplacementPlan.of({}))
    base = assess(eigen_analysis(base_case.system.A), base_case.system.labels)
    base_z0 = base_case.equilibrium.z

    rows: list[dict[str, object]] = []
    for policy in POLICIES:
        for bus in network.replacement_candidates:
            for rho in RHO_GRID:
                plan = ReplacementPlan.of({bus: rho}, q_policy=policy)
                row = {
                    "bus": bus,
                    "rho": rho,
                    "q_policy": policy,
                    "scr": metrics[bus].scr,
                    "thevenin_magnitude_pu": metrics[bus].thevenin_magnitude,
                    "rating_mva": metrics[bus].rating_mva,
                }
                try:
                    row.update(evaluate(plan, base_z0 if policy == "matched" else None))
                    row["status"] = "SUCCESS"
                    row["delta_zeta"] = row["zeta_min"] - base.zeta_min
                except InfeasibleReplacement as error:
                    row["status"] = "INFEASIBLE"
                    row["reason"] = str(error)
                rows.append(row)
    table = pd.DataFrame(rows)

    # rho* : largest feasible fraction that neither destabilizes nor degrades
    # damping materially, per bus and policy.
    limits = []
    for policy in POLICIES:
        for bus in network.replacement_candidates:
            subset = table[
                (table.bus == bus)
                & (table.q_policy == policy)
                & (table.status == "SUCCESS")
            ].sort_values("rho")
            allowed = 0.0
            for _, r in subset.iterrows():
                if bool(r.unstable) or r.delta_zeta < -MATERIAL_DEGRADATION:
                    break
                allowed = float(r.rho)
            best = subset[subset.rho == allowed] if allowed > 0 else subset.iloc[:0]
            limits.append(
                {
                    "bus": bus,
                    "q_policy": policy,
                    "rho_star": allowed,
                    "scr": metrics[bus].scr,
                    "rating_mva": metrics[bus].rating_mva,
                    "mw_replaceable": allowed * metrics[bus].rating_mva,
                    "zeta_at_rho_star": float(best.zeta_min.iloc[0])
                    if len(best)
                    else np.nan,
                    "delta_zeta_at_full": float(
                        subset[subset.rho == 1.0].delta_zeta.iloc[0]
                    )
                    if len(subset[subset.rho == 1.0])
                    else np.nan,
                    "feasible_cases": int(len(subset)),
                }
            )
    limit_table = pd.DataFrame(limits)

    matched = limit_table[limit_table.q_policy == "matched"]
    full = table[
        (table.q_policy == "matched") & (table.rho == 1.0) & (table.status == "SUCCESS")
    ]
    correlations = {}
    if len(full) > 2:
        for metric in ("scr", "thevenin_magnitude_pu"):
            rho_, p = spearmanr(full[metric], full["delta_zeta"])
            correlations[f"spearman_delta_zeta_vs_{metric}"] = {
                "rho": float(rho_),
                "p_value": float(p),
            }

    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_atlas.csv", index=False)
    limit_table.to_csv(tables / f"{EXPERIMENT}_limits.csv", index=False)

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "rho_grid": list(RHO_GRID),
            "policies": list(POLICIES),
            "material_degradation": MATERIAL_DEGRADATION,
            "base_zeta_min": base.zeta_min,
        },
    )
    manifest.finish(
        "SUCCESS",
        cases=len(table),
        infeasible=int((table.status == "INFEASIBLE").sum()),
        max_operating_point_drift_matched=float(
            table[table.q_policy == "matched"].operating_point_drift.max()
        ),
        base_zeta_min=base.zeta_min,
        base_spectral_abscissa=base.spectral_abscissa,
        rho_star_matched={int(r.bus): float(r.rho_star) for _, r in matched.iterrows()},
        correlations=correlations,
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "base: zeta_min = %+.6f at %.4f Hz,  abscissa = %+.6f"
        % (base.zeta_min, base.zeta_min_frequency_hz, base.spectral_abscissa)
    )
    print(
        "max operating-point drift under the matched policy: %.3e"
        % table[table.q_policy == "matched"].operating_point_drift.max()
    )
    print()
    for policy in POLICIES:
        print(f"--- policy: {policy}")
        view = table[table.q_policy == policy].pivot_table(
            index="bus", columns="rho", values="zeta_min"
        )
        print(view.to_string(float_format=lambda v: f"{v:8.5f}", na_rep="  INFEAS"))
        print()
    print("rho* (largest fraction with no instability and no material degradation)")
    print(
        limit_table[
            [
                "bus",
                "q_policy",
                "scr",
                "rating_mva",
                "rho_star",
                "mw_replaceable",
                "delta_zeta_at_full",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:10.4f}", na_rep="-")
    )
    print()
    print("nodal ranking vs actual damping change at full replacement (matched):")
    for name, value in correlations.items():
        print(f"  {name}: rho = {value['rho']:+.4f}  p = {value['p_value']:.4f}")
    print()
    print(f"{len(table)} cases in {time.time() - started:.1f} s")
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
