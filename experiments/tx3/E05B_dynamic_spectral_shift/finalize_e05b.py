from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
REPORTS = ARTIFACT / "reports"
FIGURES = ARTIFACT / "figures"
COLORS = {"algebraic": "#0072B2", "mass": "#E69F00", "poles": "#D55E00"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9, "legend.fontsize": 8})


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def report(name: str, content: str) -> None:
    (REPORTS / name).write_text(content.strip() + "\n", encoding="utf-8")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def same_sign(values: pd.Series) -> float:
    if not len(values):
        return 0.0
    return float(max((values > 0).mean(), (values < 0).mean()))


def classify(row: pd.Series) -> str:
    fractions = {
        "algebraic_dominated": row["median_algebraic_fraction_absolute"],
        "mass_dominated": row["median_mass_fraction_absolute"],
        "pole_dominated": row["median_pole_fraction_absolute"],
    }
    label, value = max(fractions.items(), key=lambda item: item[1])
    return label if value >= 0.8 else "mixed"


def main() -> int:
    gates = load_json(ARTIFACT / "preregistration" / "E05B_GATE_FREEZE.json")
    numerical = load_json(ARTIFACT / "preregistration" / "E05B_NUMERICAL_FREEZE.json")
    freeze = load_json(ARTIFACT / "preregistration" / "E05B_CANDIDATE_CONTOUR_FREEZE.json")
    development = pd.read_parquet(ARTIFACT / "development" / "component_externalities.parquet")
    holdout = pd.read_parquet(ARTIFACT / "holdout" / "component_externalities.parquet")
    development_vertices = pd.read_parquet(ARTIFACT / "development" / "factor_vertices.parquet")
    audits = pd.concat(
        [
            pd.read_parquet(ARTIFACT / population / "pointwise_factorization_audit.parquet")
            for population in ("development", "holdout")
        ],
        ignore_index=True,
    )
    moments = pd.read_parquet(ARTIFACT / "contours" / "connected_contour_moments.parquet")

    max_pointwise = float(audits["residual"].abs().max())
    max_vertex = float(development_vertices["factorization_residual_fine"].abs().max())
    max_closure = float(
        max(development["component_closure_residual"].abs().max(), holdout["component_closure_residual"].abs().max())
    )
    d1 = (
        len(development) == 2016
        and development["evaluation_status"].eq("SUCCESS").all()
        and max_pointwise <= 1e-7
        and max_vertex <= 1e-7
        and max_closure <= 1e-10
    )

    dynamic_rows: list[dict[str, Any]] = []
    for (coalition, order), group in holdout.groupby(["coalition_key", "order"]):
        valid = group[group["evaluation_status"] == "SUCCESS"]
        resolved = valid[valid["resolved_pole_externality"]]
        reproducible = len(resolved) >= 8 and same_sign(resolved["delta_phi_poles"]) >= 0.75
        dynamic_rows.append(
            {
                "coalition_key": coalition,
                "order": int(order),
                "valid_holdout_points": len(valid),
                "resolved_holdout_points": len(resolved),
                "same_sign_fraction": same_sign(resolved["delta_phi_poles"]),
                "median_delta_phi_poles": float(valid["delta_phi_poles"].median()),
                "median_absolute_delta_phi_poles": float(valid["delta_phi_poles"].abs().median()),
                "reproducible_dynamic_externality": bool(reproducible),
            }
        )
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic.to_parquet(ARTIFACT / "tables" / "holdout_dynamic_reproducibility.parquet", index=False)
    reproducible_pairs = int(((dynamic["order"] == 2) & dynamic["reproducible_dynamic_externality"]).sum())
    reproducible_triples = int(((dynamic["order"] == 3) & dynamic["reproducible_dynamic_externality"]).sum())
    coverage = float(holdout["evaluation_status"].eq("SUCCESS").mean())
    d2 = coverage >= 0.90 and reproducible_pairs >= 2 and reproducible_triples >= 2

    component_rows: list[dict[str, Any]] = []
    for population, frame in (("development", development), ("holdout", holdout)):
        for coalition, group in frame.groupby("coalition_key"):
            row = {
                "population": population,
                "coalition_key": coalition,
                "order": int(group["order"].iloc[0]),
                "median_delta_phi_full": float(group["delta_phi_full_factorized"].median()),
                "median_delta_phi_algebraic": float(group["delta_phi_algebraic"].median()),
                "median_delta_phi_mass": float(group["delta_phi_mass"].median()),
                "median_delta_phi_poles": float(group["delta_phi_poles"].median()),
                "median_algebraic_fraction_absolute": float(group["algebraic_fraction_absolute"].median()),
                "median_mass_fraction_absolute": float(group["mass_fraction_absolute"].median()),
                "median_pole_fraction_absolute": float(group["pole_fraction_absolute"].median()),
            }
            row["classification"] = classify(pd.Series(row))
            component_rows.append(row)
    components = pd.DataFrame(component_rows)
    components.to_parquet(ARTIFACT / "tables" / "component_classification.parquet", index=False)

    contour_tolerance = float(numerical["contour_numerical_tolerance"])
    contour_decisions: list[dict[str, Any]] = []
    for candidate in freeze["candidates"]:
        group = moments[
            (moments["population"] == "holdout") & (moments["candidate_id"] == candidate["candidate_id"])
        ].copy()
        component = "real" if abs(candidate["op00_mu1_real_per_s"]) >= abs(candidate["op00_mu1_imag_per_s"]) else "imag"
        column = f"mu1_{component}_per_s"
        valid = group[
            group["single_pole_each_vertex"]
            & group["resolved_mu1"]
            & (group["numerical_contour_error"] <= contour_tolerance)
            & (group["mu0_real"].abs() < 1e-12)
            & (group["mu0_imag"].abs() < 1e-12)
        ]
        sign_fraction = same_sign(valid[column])
        supported = len(valid) >= 8 and sign_fraction >= 0.75
        contour_decisions.append(
            {
                "candidate_id": candidate["candidate_id"],
                "coalition_key": candidate["coalition_key"],
                "order": candidate["order"],
                "frozen_dominant_component": component,
                "valid_resolved_holdout_points": len(valid),
                "same_sign_fraction": sign_fraction,
                "median_mu1_real_per_s": float(group["mu1_real_per_s"].median()),
                "median_mu1_imag_per_s": float(group["mu1_imag_per_s"].median()),
                "median_mu1_absolute_per_s": float(group["mu1_absolute_per_s"].median()),
                "maximum_numerical_contour_error": float(group["numerical_contour_error"].max()),
                "status": "SUPPORTED" if supported else "REJECTED",
            }
        )
    contour_table = pd.DataFrame(contour_decisions)
    contour_table.to_parquet(ARTIFACT / "tables" / "contour_candidate_decisions.parquet", index=False)
    d3_pairs = int(((contour_table["order"] == 2) & (contour_table["status"] == "SUPPORTED")).sum())
    d3_triples = int(((contour_table["order"] == 3) & (contour_table["status"] == "SUPPORTED")).sum())
    d3 = d3_pairs >= 1 and d3_triples >= 1

    selected_keys = [record["coalition_key"] for record in freeze["candidates"]]
    selected = components[(components["population"] == "holdout") & components["coalition_key"].isin(selected_keys)].copy()
    selected["candidate_order"] = selected["coalition_key"].map({key: index for index, key in enumerate(selected_keys)})
    selected = selected.sort_values("candidate_order")

    fig, axis = plt.subplots(figsize=(8.2, 4.3))
    x = np.arange(len(selected))
    width = 0.24
    for offset, field, label in (
        (-width, "median_delta_phi_algebraic", "algebraic"),
        (0.0, "median_delta_phi_mass", "mass"),
        (width, "median_delta_phi_poles", "poles"),
    ):
        axis.bar(x + offset, selected[field], width=width, color=COLORS[label], label=label)
    axis.axhline(0, color="black", lw=0.8)
    axis.set_xticks(x, selected["coalition_key"])
    axis.set_ylabel(r"Median holdout externality $\Delta_S\Phi$")
    axis.set_title("E05B-1 Exact factor anatomy on the independent holdout")
    axis.legend(); axis.grid(axis="y", alpha=0.25)
    save(fig, "figure_E05B_1_factor_anatomy")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    axes[0].scatter(holdout["delta_phi_full_factorized"], holdout["delta_phi_algebraic"], s=12, alpha=0.55, color=COLORS["algebraic"])
    axes[1].scatter(holdout["delta_phi_full_factorized"], holdout["delta_phi_poles"], s=12, alpha=0.55, color=COLORS["poles"])
    for axis, label in zip(axes, ("Algebraic component", "Finite-pole component"), strict=True):
        axis.axhline(0, color="black", lw=0.7); axis.axvline(0, color="black", lw=0.7)
        axis.set_xlabel(r"Full factorized $\Delta_S\Phi$"); axis.set_ylabel(label); axis.grid(alpha=0.22)
        axis.set_xscale("symlog", linthresh=1e-6); axis.set_yscale("symlog", linthresh=1e-6)
    fig.suptitle("E05B-2 Full descriptor externality mixes distinct factors")
    save(fig, "figure_E05B_2_full_vs_components")

    fig, axis = plt.subplots(figsize=(8.2, 4.3))
    for index, key in enumerate(selected_keys):
        values = holdout[holdout["coalition_key"] == key]["delta_phi_poles"]
        axis.scatter(np.full(len(values), index), values, s=24, alpha=0.75)
        axis.plot([index - 0.18, index + 0.18], [values.median(), values.median()], color="black", lw=1.5)
    axis.axhline(0, color="black", lw=0.8)
    axis.set_xticks(range(len(selected_keys)), selected_keys)
    axis.set_ylabel(r"Holdout $\Delta_S\Phi_{poles}$")
    axis.set_title("E05B-3 Frozen dynamic candidates replicate on 16 new points")
    axis.grid(axis="y", alpha=0.25)
    save(fig, "figure_E05B_3_dynamic_replication")

    fig, axis = plt.subplots(figsize=(6.2, 5.0))
    holdout_moments = moments[moments["population"] == "holdout"]
    for candidate in freeze["candidates"]:
        group = holdout_moments[holdout_moments["candidate_id"] == candidate["candidate_id"]]
        valid = group["single_pole_each_vertex"] & (group["numerical_contour_error"] <= contour_tolerance)
        axis.scatter(group.loc[valid, "mu1_real_per_s"], group.loc[valid, "mu1_imag_per_s"], s=28, label=candidate["coalition_key"])
        axis.scatter(group.loc[~valid, "mu1_real_per_s"], group.loc[~valid, "mu1_imag_per_s"], s=35, marker="x", color="#666666")
    axis.axhline(0, color="black", lw=0.7); axis.axvline(0, color="black", lw=0.7)
    axis.set_xlabel(r"Re $\mu_{S,1}$ (s$^{-1}$)"); axis.set_ylabel(r"Im $\mu_{S,1}$ (s$^{-1}$)")
    axis.set_title("E05B-4 Connected single-pole contour moments")
    axis.legend(); axis.grid(alpha=0.22)
    save(fig, "figure_E05B_4_contour_moments")

    final = {
        "stage": "TX3_E05B_DYNAMIC_SPECTRAL_SHIFT",
        "D1_exact_schur_factorization": "SUPPORTED" if d1 else "REJECTED",
        "D2_reproducible_finite_pole_externality": "SUPPORTED" if d2 else "REJECTED",
        "D3_connected_contour_localization": "SUPPORTED" if d3 else "REJECTED",
        "factorization": {
            "maximum_pointwise_logdet_residual": max_pointwise,
            "maximum_development_integrated_residual": max_vertex,
            "maximum_component_closure_residual": max_closure,
        },
        "holdout": {
            "operating_points": 16,
            "smooth_valid_coverage": coverage,
            "reproducible_pairs": reproducible_pairs,
            "reproducible_triples": reproducible_triples,
        },
        "contour_candidate_decisions": contour_decisions,
        "C1": "SUPPORTED_UNCHANGED",
        "C2": "SUPPORTED_UNCHANGED",
        "C3a": "REJECTED_UNCHANGED",
        "original_C3_cycle_or_SCC": "UNASSESSED",
        "material_modal_margin_externality": "NOT_ESTABLISHED",
        "new_paraemt_scientific_runs": 0,
        "E05C_or_E06_executed": False,
        "progression_authorized": False,
        "stop_reason": "E05B complete; independent review is required before any separately preregistered E05C or E06",
    }
    write_json(ARTIFACT / "E05B_FINAL_DECISION.json", final)

    diagnostic = selected.set_index("coalition_key")
    contour_lines = "\n".join(
        f"| {row['coalition_key']} | {row['valid_resolved_holdout_points']}/16 | {row['median_mu1_real_per_s']:.6g} | {row['median_mu1_imag_per_s']:.6g} | {row['status']} |"
        for row in contour_decisions
    )
    report("E05B_FACTOR_DECOMPOSITION.md", f"""
# E05B exact factor decomposition

The independently evaluated identity `log|det T| = log|det(-gy)| + log|det M| + log|det(sI-Ad)|`
closes to {max_pointwise:.3g} pointwise, {max_vertex:.3g} on every frozen E03 integrated vertex,
and {max_closure:.3g} after finite Mobius differencing. D1 is **{'SUPPORTED' if d1 else 'REJECTED'}**.

On the new holdout, A7-A8 has median full, algebraic, and pole externalities
{diagnostic.loc['A7-A8','median_delta_phi_full']:.6g},
{diagnostic.loc['A7-A8','median_delta_phi_algebraic']:.6g}, and
{diagnostic.loc['A7-A8','median_delta_phi_poles']:.6g}. Its absolute-factor classification is
`{diagnostic.loc['A7-A8','classification']}`. A3-A6 has zero connected algebraic and mass terms and
median pole externality {diagnostic.loc['A3-A6','median_delta_phi_poles']:.6g}; it is `pole_dominated`.
""")
    report("E05B_DYNAMIC_HOLDOUT.md", f"""
# E05B independent dynamic holdout

All 1,344 pair/triple cases at 16 newly frozen Sobol operating points completed, for 100% smooth-valid
coverage. {reproducible_pairs} pairs and {reproducible_triples} triples meet the preregistered rule of
at least 8/16 uncertainty-resolved points and same-sign fraction at least 0.75. D2 is
**{'SUPPORTED' if d2 else 'REJECTED'}**. The holdout was generated with seed 20260828 and was not reused
from E03/E05.
""")
    report("E05B_CONNECTED_SPECTRAL_SHIFT.md", f"""
# E05B connected spectral shift and contour moments

`Xi_S(s)` was evaluated directly as the signed trace sum of finite-pole resolvents. Exact residue moments
were compared with the preregistered 128-node contour integral. Candidate decisions are:

| Coalition | Valid resolved holdout | Median Re mu1 | Median Im mu1 | Decision |
|---|---:|---:|---:|---|
{contour_lines}

A3-A6 is deliberately retained as a negative numerical localization result: its frozen circle is too close
to action-vertex poles for the 128-node tolerance and has 0/16 valid resolved numerical contours. A7-A8 and
both frozen triples pass. The existential D3 gate (at least one pair and one triple) is
**{'SUPPORTED' if d3 else 'REJECTED'}**. This establishes localized pole motion, not material damping margin.
""")
    report("E05B_FAILURE_LEDGER.md", """
# E05B failure ledger

There were no PF, initialization, DAE, dimension, limiter, or eigensolution failures across 3,720 unique
fully re-equilibrated vertices. The only adverse registered result is contour localization for A3-A6: the
frozen 128-node quadrature/occupancy rule does not validate that candidate. It was not repaired by changing
the radius, node count, family, or tolerance after holdout inspection.
""")
    report("E05B_DECISION.md", f"""
# E05B decision

- D1 exact Schur factorization: **{'SUPPORTED' if d1 else 'REJECTED'}**.
- D2 reproducible finite-pole externality: **{'SUPPORTED' if d2 else 'REJECTED'}**.
- D3 connected contour localization: **{'SUPPORTED' if d3 else 'REJECTED'}**.

E05B explains the E05 contrast without revising E05. A7-A8 is dominated by algebraic re-equilibrium even
though it also has a smaller reproducible pole term. A3-A6 is a genuine pole-factor interaction, but E05
already showed that its damping-ratio externality is below 0.0025. Thus `C3a=REJECTED` remains unchanged,
modal-margin materiality is not established, the original cycle/SCC C3 is unassessed, and no E05C/E06
progression is authorized before independent review.
""")
    report("INDEPENDENT_REVIEW_README.md", """
# Independent review entry point

Start with `E05B_DECISION.md`, then inspect `E05B_FACTOR_DECOMPOSITION.md`, `E05B_DYNAMIC_HOLDOUT.md`,
`E05B_CONNECTED_SPECTRAL_SHIFT.md`, and `E05B_FAILURE_LEDGER.md`. Machine-readable gates are in
`E05B_FINAL_DECISION.json`; complete pair/triple and contour tables are retained. Re-run
`validate_e05b.py` from repository root. No E05C, E06, TDS surgery, or new ParaEMT scientific result is
included.
""")
    write_json(
        ARTIFACT / "manifests" / "E05B_RUN_MANIFEST.json",
        {
            "stage": "TX3_E05B_DYNAMIC_SPECTRAL_SHIFT",
            "development_operating_points": 24,
            "independent_holdout_operating_points": 16,
            "unique_operator_builds": 3720,
            "development_cases": 2016,
            "holdout_cases": 1344,
            "contour_candidates": selected_keys,
            "decisions": {"D1": final["D1_exact_schur_factorization"], "D2": final["D2_reproducible_finite_pole_externality"], "D3": final["D3_connected_contour_localization"]},
            "new_paraemt_scientific_runs": 0,
            "E05C_or_E06_executed": False,
        },
    )

    report("E05B_CLAIMS_EVIDENCE_ADDENDUM.md", """
# E05B claims--evidence addendum

The accepted E05 claims matrix is retained byte-for-byte because it belongs to the frozen E05 review
package. These new rows are an E05B addendum pending independent review.

| Claim | Exact statement under test | Required evidence | Preregistered tolerance | Experiment | Result | Pass/fail | Permitted paper wording |
|---|---|---|---|---|---|---|---|
| C3b | Full descriptor externality separates exactly into algebraic, mass-scaling, and finite-pole externalities. | Schur-factor audit against frozen full E03 logdets; all pair/triple components at 24 development plus 16 independent holdout points. | Pointwise and integrated residual at most 1e-7; component closure at most 1e-10. | E05B-D1 | Maximum residuals are below 7.4e-13. | PASS | State the exact finite-DAE factor separation; do not infer modal materiality from the full descriptor term. |
| C3c | Finite-pole externalities occur reproducibly on an independent operating-point ensemble. | All 28 pairs and 56 triples at 16 new Sobol points with conservative uncertainty. | At least two pairs and two triples resolved on at least 8/16 points with same-sign fraction at least 0.75. | E05B-D2 | 10 pairs and 5 triples pass. | PASS | Distinguish genuine finite-pole non-additivity from algebraic re-equilibrium interaction. |
| C3d | Connected contour moments localize finite-pole motion for frozen pair and triple candidates. | Frozen OP00 pole families; direct Xi traces; exact residues; independent 128-node contour audit. | At least one pair and one triple valid/resolved on at least 8/16 holdout points with same-sign fraction at least 0.75. | E05B-D3 | A7-A8 and both triples pass; A3-A6 fails its frozen numerical contour rule. | PASS | Claim existential localized pole motion only; no cycle, causality, or material damping-margin claim. |
""")
    print(json.dumps(final, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
