"""E12 - PHASE C. Replacement compatibility census (R2).

Question
    Do singles and pairs allow the outcome of a larger replacement portfolio to
    be inferred, and where they do not, is the failure a genuine high-order
    interaction or the accumulation of lower-order ones crossing a threshold?

Independent variables
    every subset of the nine candidate buses up to size three, at rho = 1.

Dependent variables
    instability, minimum damping in 0.1-5 Hz, damping change against base,
    replaced MW, and the Moebius interaction orders of both functionals.

Controls
    every subset is evaluated by the same pipeline; the mechanism classifier is
    automated and cannot be overridden by inspection.

Status of the result
    NUMERICAL OBSERVATION at one operating point, one converter tuning and
    rho = 1.

Usage
    python experiments/E12_compatibility_census.py [--max-size 3]
"""

from __future__ import annotations

import argparse
import time
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network

EXPERIMENT = "E12_compatibility_census"
MATERIAL_DEGRADATION = 0.005
POLICY = "matched"


def classify(
    values: dict[frozenset[str], float],
    identity: dict[frozenset[str], str],
    members: tuple[str, ...],
) -> tuple[str, dict[str, float]]:
    """Genuine high-order effect, or lower-order accumulation crossing zero?

    The functional is the spectral abscissa over non-reference modes, because its
    SIGN is the criterion. Using minimum damping in a band instead would compare
    a decomposition of one quantity against a verdict taken from another.

    The abscissa is a maximum over modes and therefore not smooth: if the mode
    attaining it changes across the subset lattice, the Moebius terms mix
    different physical modes and no interaction order can be claimed. That case
    is reported as C_MODE_SWITCHING rather than forced into A or B, which is what
    the preregistration reserves class C for.
    """

    decomposition = decompose(lambda s: values[frozenset(s)], members)
    order = len(members)
    truncations = [decomposition.truncated(members, k) for k in range(order + 1)]
    exact = truncations[-1]
    evidence = {f"order_{k}": float(value) for k, value in enumerate(truncations)}
    evidence["irreducible_top"] = float(decomposition.terms[frozenset(members)])
    families = {
        identity[frozenset(sub)]
        for size in range(order + 1)
        for sub in combinations(members, size)
    }
    evidence["distinct_critical_families"] = float(len(families))
    if exact <= 0.0:
        return "SAFE", evidence
    if len(families) > 1:
        return "C_MODE_SWITCHING", evidence
    lower_crosses = any(value > 0.0 for value in truncations[:-1])
    if not lower_crosses:
        return "A_GENUINE_HIGH_ORDER", evidence
    return "B_LOWER_ORDER_ACCUMULATION", evidence


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE C compatibility census")
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--max-size", type=int, default=9)
    args = parser.parse_args()
    started = time.time()

    network = load_network()
    candidates = network.replacement_candidates

    records: dict[tuple[int, ...], dict[str, object]] = {}
    failures: list[dict[str, object]] = []
    for size in range(args.max_size + 1):
        for members in combinations(candidates, size):
            plan = ReplacementPlan.of({b: 1.0 for b in members}, q_policy=POLICY)
            try:
                case = solve_case(plan)
                spectrum = eigen_analysis(case.system.A)
                verdict = assess(spectrum, case.system.labels)
                dynamic = [m for m in spectrum.modes if abs(m.value) > 1e-3]
                critical = max(dynamic, key=lambda m: m.real)
                label = case.system.labels[int(np.argmax(critical.participation))]
                family = "gfl" if "gfl" in label else "sg" if "sg" in label else "other"
                records[members] = {
                    "critical_dominant_state": label,
                    "critical_family": family,
                    "members": "+".join(map(str, members)) or "BASE",
                    "size": size,
                    "replaced_mw": case.replaced_mw,
                    "n_states": case.n_states,
                    "spectral_abscissa": verdict.spectral_abscissa,
                    "zeta_min": verdict.zeta_min,
                    "zeta_min_frequency_hz": verdict.zeta_min_frequency_hz,
                    "critical_frequency_hz": verdict.critical_frequency_hz,
                    "critical_condition": verdict.critical_condition,
                    "rhp_modes": verdict.rhp_modes,
                    "unstable": verdict.unstable,
                    "gz_condition": case.gz_condition,
                    "equilibrium_norm_g": case.equilibrium.norm_g,
                }
            except InfeasibleReplacement as error:
                failures.append({"members": members, "reason": str(error)})

    base = records[()]
    base_zeta = float(base["zeta_min"])
    threshold = base_zeta - MATERIAL_DEGRADATION

    zeta = {frozenset(map(str, k)): float(v["zeta_min"]) for k, v in records.items()}
    for record in records.values():
        record["delta_zeta"] = float(record["zeta_min"]) - base_zeta
        record["safe"] = (
            not record["unstable"] and float(record["zeta_min"]) >= threshold
        )

    # Non-composability analysis over every portfolio whose proper subsets are safe.
    abscissa = {
        frozenset(map(str, k)): float(v["spectral_abscissa"])
        for k, v in records.items()
    }
    identity = {
        frozenset(map(str, k)): str(v["critical_family"]) for k, v in records.items()
    }
    findings: list[dict[str, object]] = []
    degradations: list[dict[str, object]] = []
    for members, record in records.items():
        if len(members) < 2:
            continue
        subsets_stable = all(
            not records[sub]["unstable"]
            for size in range(1, len(members))
            for sub in combinations(members, size)
        )
        subsets_safe = all(
            records[sub]["safe"]
            for size in range(1, len(members))
            for sub in combinations(members, size)
        )
        record["all_proper_subsets_safe"] = subsets_safe
        record["all_proper_subsets_stable"] = subsets_stable
        names = tuple(map(str, members))
        if subsets_stable and record["unstable"]:
            mechanism, evidence = classify(abscissa, identity, names)
            findings.append(
                {
                    "members": record["members"],
                    "size": len(members),
                    "spectral_abscissa": record["spectral_abscissa"],
                    "zeta_min": record["zeta_min"],
                    "critical_frequency_hz": record["critical_frequency_hz"],
                    "critical_family": record["critical_family"],
                    "replaced_mw": record["replaced_mw"],
                    "mechanism": mechanism,
                    **evidence,
                }
            )
        elif subsets_safe and not record["safe"]:
            degradations.append(
                {
                    "members": record["members"],
                    "size": len(members),
                    "zeta_min": record["zeta_min"],
                    "delta_zeta": record["delta_zeta"],
                    "unstable": record["unstable"],
                    "replaced_mw": record["replaced_mw"],
                }
            )

    table = pd.DataFrame(list(records.values()))
    # Additivity check: how well does the pairwise model predict triples?
    prediction_rows = []
    for members, record in records.items():
        if len(members) < 2:
            continue
        names = tuple(map(str, members))
        decomposition = decompose(lambda s: zeta[frozenset(s)], names)
        row = {
            "members": record["members"],
            "size": len(members),
            "exact": record["zeta_min"],
            "order_1": decomposition.truncated(names, 1),
            "order_2": decomposition.truncated(names, 2),
            "irreducible_top": decomposition.terms[frozenset(names)],
            "additive_error": decomposition.truncated(names, 1) - record["zeta_min"],
            "pairwise_error": decomposition.truncated(names, 2) - record["zeta_min"],
        }
        if len(members) >= 3:
            row["order_3"] = decomposition.truncated(names, 3)
            row["triplewise_error"] = row["order_3"] - record["zeta_min"]
        prediction_rows.append(row)
    predictions = pd.DataFrame(prediction_rows)

    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_portfolios.csv", index=False)
    predictions.to_csv(tables / f"{EXPERIMENT}_triple_predictions.csv", index=False)
    if findings:
        pd.DataFrame(findings).to_csv(
            tables / f"{EXPERIMENT}_non_composable.csv", index=False
        )
    if degradations:
        pd.DataFrame(degradations).to_csv(
            tables / f"{EXPERIMENT}_degradation_only.csv", index=False
        )

    by_size = table.groupby("size").agg(
        cases=("members", "count"),
        unstable=("unstable", "sum"),
        unsafe=("safe", lambda s: int((~s).sum())),
        zeta_min_worst=("zeta_min", "min"),
        zeta_min_best=("zeta_min", "max"),
    )

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "policy": POLICY,
            "rho": 1.0,
            "max_size": args.max_size,
            "material_degradation": MATERIAL_DEGRADATION,
            "candidates": list(candidates),
        },
    )
    manifest.finish(
        "SUCCESS",
        cases=len(records),
        infeasible=len(failures),
        base_zeta_min=base_zeta,
        safety_threshold=threshold,
        unstable_portfolios=int(table.unstable.sum()),
        unsafe_portfolios=int((~table.safe).sum()),
        non_composable=findings,
        degradation_only=degradations,
        minimum_unstable_order=(
            int(table[table.unstable]["size"].min())
            if bool(table.unstable.any())
            else None
        ),
        mechanism_counts=(
            pd.DataFrame(findings).mechanism.value_counts().to_dict()
            if findings
            else {}
        ),
        worst_triple=(
            predictions.sort_values("exact").iloc[0].to_dict()
            if len(predictions)
            else None
        ),
        max_additive_prediction_error=(
            float(predictions.additive_error.abs().max()) if len(predictions) else None
        ),
        max_pairwise_prediction_error=(
            float(predictions.pairwise_error.abs().max()) if len(predictions) else None
        ),
        truncation_error_by_size=(
            predictions.groupby("size")
            .agg(
                additive_max=("additive_error", lambda s: float(s.abs().max())),
                pairwise_max=("pairwise_error", lambda s: float(s.abs().max())),
            )
            .to_dict("index")
            if len(predictions)
            else {}
        ),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("base zeta_min = %+.6f   safety threshold = %+.6f" % (base_zeta, threshold))
    print()
    print(by_size.to_string())
    print()
    if findings:
        found = pd.DataFrame(findings)
        print("MINIMAL INCOMPATIBLE CORES: every proper subset stable, portfolio not")
        print("  count = %d   smallest order = %d" % (len(found), found["size"].min()))
        print("  mechanism: %s" % found.mechanism.value_counts().to_dict())
        print(
            found[found["size"] == found["size"].min()][
                [
                    "members",
                    "spectral_abscissa",
                    "critical_frequency_hz",
                    "critical_family",
                    "replaced_mw",
                    "mechanism",
                    "distinct_critical_families",
                ]
            ].to_string(index=False)
        )
    else:
        print("MINIMAL INCOMPATIBLE CORES: none found up to size %d" % args.max_size)
    print()
    print(
        "degradation-only portfolios (stable but below the damping threshold): %d"
        % len(degradations)
    )
    print()
    print("worst portfolios by zeta_min:")
    print(
        table.nsmallest(8, "zeta_min")[
            [
                "members",
                "size",
                "zeta_min",
                "delta_zeta",
                "replaced_mw",
                "safe",
                "unstable",
            ]
        ].to_string(index=False)
    )
    print()
    if len(predictions):
        print("truncation error of a lower-order model, by portfolio size")
        summary = predictions.groupby("size").agg(
            cases=("members", "count"),
            additive_max=("additive_error", lambda s: s.abs().max()),
            pairwise_max=("pairwise_error", lambda s: s.abs().max()),
            pairwise_median=("pairwise_error", lambda s: s.abs().median()),
            irreducible_top_max=("irreducible_top", lambda s: s.abs().max()),
        )
        print(summary.to_string(float_format=lambda v: f"{v:10.6f}"))
        print()
        print("worst pairwise-model errors overall:")
        print(
            predictions.reindex(
                predictions.pairwise_error.abs().sort_values(ascending=False).index
            )
            .head(8)[
                [
                    "members",
                    "size",
                    "exact",
                    "order_1",
                    "order_2",
                    "irreducible_top",
                    "pairwise_error",
                ]
            ]
            .to_string(index=False)
        )
    print()
    print(
        "%d cases, %d infeasible, %.1f s"
        % (len(records), len(failures), time.time() - started)
    )
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
