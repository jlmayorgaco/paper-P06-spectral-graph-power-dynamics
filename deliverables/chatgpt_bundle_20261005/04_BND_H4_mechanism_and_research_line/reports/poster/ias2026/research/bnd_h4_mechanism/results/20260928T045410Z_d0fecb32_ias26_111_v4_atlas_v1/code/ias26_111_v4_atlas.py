"""Read-only post-hoc extraction for IAS26-111 (frozen IAS26-060 V4 outputs).

This script never runs power flow, DAE, or eigenvalue solves and never writes
to the frozen source run. All outputs are scoped to this immutable ticket run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
import numpy as np
import pandas as pd


RUN_ID = "20260928T045410Z_d0fecb32_ias26_111_v4_atlas_v1"
HERE = Path(__file__).resolve()
RUN_ROOT = HERE.parents[1]
REPO = HERE.parents[8]
CONTRACT = REPO / "reports/poster/ias2026/research/bnd_h4_mechanism"
SOURCE_RUN_ID = "20260927T144629Z_d0fecb32_ias26_060_operating_v1"
SOURCE_RUN = CONTRACT / "results" / SOURCE_RUN_ID
MC_CASES = SOURCE_RUN / "derived/MC_CASES.csv"
SCENARIO_METRICS = SOURCE_RUN / "derived/SCENARIO_METRICS.csv"
CONFIG = CONTRACT / "configs/IAS26-060_MC_OPERATING_V1.json"
SCENARIOS = CONTRACT / "inputs/IAS26-060_SCENARIOS_V1.csv"
SOURCE_MEMO = CONTRACT / "docs/IAS26-060_V4_BLOCKER_ATLAS_DESCRIPTIVE.md"
TAU_DEC = 1e-8
BUSES = (30, 33, 35, 37)
EXPECTED = {"none": 412, "h4_only": 190, "alternative": 365, "valid": 967}
MASKS = tuple(range(16))
COLORS = {
    "NONE": "#7A8793",
    "H4_ONLY": "#D98E04",
    "{30,33,35}": "#168C83",
    "{30,33,35} | {30,33,37}": "#3976A8",
    "{30,33,35} | {30,33,37} | {33,35,37}": "#8967A5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def classify(value: float) -> str:
    if not math.isfinite(value):
        return "INVALID"
    if value < -TAU_DEC:
        return "STABLE"
    if value > TAU_DEC:
        return "UNSTABLE"
    return "INDETERMINATE"


def mask_name(mask: int) -> str:
    members = [str(bus) for bit, bus in enumerate(BUSES) if mask & (1 << bit)]
    return "{" + ",".join(members) + "}" if members else "{}"


def mask_list(raw: object) -> tuple[int, ...]:
    value = str(raw)
    if value in {"NONE", "nan", "", "None"}:
        return ()
    return tuple(sorted(int(part) for part in value.split(";") if part))


def exact_blocker_label(blockers: tuple[int, ...]) -> str:
    if not blockers:
        return "NONE"
    if blockers == (15,):
        return "H4_ONLY"
    return " | ".join(mask_name(mask) for mask in blockers)


def source_hashes() -> dict[str, str]:
    inputs = {
        "frozen_config": CONFIG,
        "frozen_scenario_manifest": SCENARIOS,
        "frozen_MC_CASES": MC_CASES,
        "frozen_SCENARIO_METRICS": SCENARIO_METRICS,
        "prior_read_only_atlas_memo": SOURCE_MEMO,
    }
    return {key: sha256(path) for key, path in inputs.items()}


def load_and_analyze() -> dict[str, object]:
    cases_raw = pd.read_csv(MC_CASES)
    metrics_raw = pd.read_csv(SCENARIO_METRICS)
    manifest_raw = pd.read_csv(SCENARIOS)
    frozen_config = json.loads(CONFIG.read_text(encoding="utf-8"))
    cases_raw["source_row_line"] = np.arange(len(cases_raw), dtype=int) + 2
    metrics_raw["source_row_line"] = np.arange(len(metrics_raw), dtype=int) + 2

    cases = cases_raw[
        (cases_raw["treatment"] == "ORIGINAL")
        & np.isclose(cases_raw["g"].astype(float), 0.03625, rtol=0.0, atol=1e-12)
    ].copy()
    cases["portfolio_mask"] = cases["portfolio_mask"].astype(int)
    metrics = metrics_raw.copy()
    metrics["scenario_valid"] = metrics["scenario_valid"].astype(str).str.lower().eq("true")
    metrics["outcome_dependent_filtering_applied"] = (
        metrics["outcome_dependent_filtering_applied"].astype(str).str.lower().eq("true")
    )
    cases["physical_feasible"] = cases["physical_feasible"].astype(str).str.lower().eq("true")

    structural_errors: list[str] = []
    if len(cases) != 16_000:
        structural_errors.append(f"expected 16000 ORIGINAL cases at g=0.03625; found {len(cases)}")
    if cases["scenario_id"].nunique() != 1000:
        structural_errors.append(f"expected 1000 case IDs; found {cases['scenario_id'].nunique()}")
    if set(metrics["scenario_id"]) != set(cases["scenario_id"]):
        structural_errors.append("case and scenario-metric ID sets differ")
    if len(manifest_raw) != 1000 or manifest_raw["scenario_id"].duplicated().any():
        structural_errors.append("frozen manifest must contain exactly 1000 unique scenario IDs")
    if set(manifest_raw["scenario_id"]) != set(metrics["scenario_id"]):
        structural_errors.append("evaluated scenario IDs differ from the frozen manifest")
    per_id_masks = cases.groupby("scenario_id")["portfolio_mask"].agg(list)
    if not per_id_masks.map(lambda values: sorted(values) == list(MASKS)).all():
        structural_errors.append("at least one ID does not contain each exact V4 mask once")
    if metrics["scenario_id"].duplicated().any() or len(metrics) != 1000:
        structural_errors.append("scenario metrics must contain exactly 1000 unique IDs")

    valid_ids = set(metrics.loc[metrics["scenario_valid"], "scenario_id"])
    valid_case_rows = cases[cases["scenario_id"].isin(valid_ids)]
    spectra_complete = valid_case_rows["alpha_perp_complete_spectrum"].astype(str).str.lower().eq("true")
    alpha_finite = np.isfinite(valid_case_rows["alpha_perp"].astype(float))
    valid_cases_pass = valid_case_rows["status"].astype(str).eq("PASS")
    if not bool(spectra_complete.all()):
        structural_errors.append("at least one valid original case lacks a complete transverse spectrum")
    if not bool(alpha_finite.all()):
        structural_errors.append("at least one valid original alpha_perp value is nonfinite")
    if not bool(valid_cases_pass.all()):
        structural_errors.append("at least one frozen-valid original case is not marked PASS")
    metric_validity = metrics.set_index("scenario_id")["scenario_valid"]
    physical_validity = cases.groupby("scenario_id")["physical_feasible"].all()
    physical_mismatch_ids = sorted(
        sid for sid in metric_validity.index
        if bool(metric_validity.loc[sid]) != bool(physical_validity.loc[sid])
    )
    if physical_mismatch_ids:
        structural_errors.append(f"scenario_valid and physical_feasible disagree for {len(physical_mismatch_ids)} IDs")
    outcome_filter_true_count = int(metrics["outcome_dependent_filtering_applied"].sum())
    if outcome_filter_true_count:
        structural_errors.append("frozen metrics record outcome-dependent filtering")
    frozen_tau = float(frozen_config["numerical_policy"]["eigenvalues"]["tau_dec_s-1"])
    frozen_g = float(frozen_config["system"]["original_treatment_g"])
    frozen_buses = tuple(int(bus) for bus in frozen_config["system"]["h4_buses_in_frozen_order"])
    if frozen_tau != TAU_DEC:
        structural_errors.append(f"script tau_dec={TAU_DEC} differs from frozen config tau_dec={frozen_tau}")
    if frozen_g != 0.03625:
        structural_errors.append(f"script original g=0.03625 differs from frozen config g={frozen_g}")
    if frozen_buses != BUSES:
        structural_errors.append(f"script V4 port order {BUSES} differs from frozen config {frozen_buses}")
    metric_by_id = metrics.set_index("scenario_id", drop=False)
    case_groups = {sid: group.sort_values("portfolio_mask") for sid, group in cases.groupby("scenario_id", sort=False)}

    records: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    portfolio_class_counts = {mask: Counter() for mask in MASKS}
    blocker_membership = Counter()
    exact_set_counts = Counter()
    kappa_counts = Counter()
    cardinality_any_unstable = Counter()
    cardinality_any_indeterminate = Counter()
    triple_membership = Counter()
    scenario_results: dict[str, dict[str, object]] = {}
    metric_mask_mismatches: list[str] = []
    metric_xy_mismatches: list[str] = []
    base_unstable_ids: list[str] = []

    for sid, group in case_groups.items():
        row = metric_by_id.loc[sid]
        alpha = {int(r.portfolio_mask): float(r.alpha_perp) for r in group.itertuples()}
        source_lines = ";".join(str(int(x)) for x in group["source_row_line"].tolist())
        is_valid = sid in valid_ids
        raw_failures = sorted({str(v) for v in group["failure_code"].dropna() if str(v).strip()})
        nonpass = sorted({str(v) for v in group["status"].dropna() if str(v) != "PASS"})
        failure_codes = ";".join(raw_failures or nonpass)

        # Retain every raw portfolio spectrum column for all 1000 IDs. Classify
        # only the 967 physically valid scenarios admitted by the frozen gate.
        result: dict[str, object] = {
            "scenario_id": sid,
            "scenario_index_zero_based": int(row["scenario_index_zero_based"]),
            "source_run_id": SOURCE_RUN_ID,
            "scenario_valid": bool(is_valid),
            "atlas_status": "INCLUDED_VALID" if is_valid else "EXCLUDED_PHYSICAL_FEASIBILITY",
            "failure_codes": failure_codes,
            "load_scale": float(row["load_scale"]),
            "epsilon_l2": float(row["epsilon_l2"]),
            "saturation_count": int(row["saturation_count"]),
            "P_L_total_MW": float(row["P_L_total_MW"]),
            "Q_L_total_Mvar": float(row["Q_L_total_Mvar"]),
            "delta_P_G_requested_total_MW": float(row["delta_P_G_requested_total_MW"]),
            "redispatch_saturated_buses_json": row["redispatch_saturated_buses_json"],
            "outcome_dependent_filtering_applied": bool(row["outcome_dependent_filtering_applied"]),
            "x_s_max_proper_alpha_s-1": math.nan,
            "y_s_H4_alpha_s-1": math.nan,
            "m_comp_min_y_minus_x_s-1": math.nan,
            "minimal_blocker_masks": "NOT_CLASSIFIED_INVALID" if not is_valid else "",
            "minimal_blocker_ids": "NOT_CLASSIFIED_INVALID" if not is_valid else "",
            "exact_blocker_set_category": "INVALID_PHYSICAL_SCENARIO" if not is_valid else "",
            "n_minimal_blockers": math.nan,
            "kappa_minimal_blocker_order": math.nan,
            "h4_unstable": "NOT_CLASSIFIED_INVALID" if not is_valid else "",
            "metric_reconciliation": "NOT_IN_VALID_ATLAS" if not is_valid else "",
        }
        for mask in MASKS:
            result[f"alpha_perp_mask_{mask:02d}_s-1"] = alpha[mask]

        blockers: tuple[int, ...] = ()
        x_value = y_value = margin = math.nan
        if is_valid:
            states = {mask: classify(value) for mask, value in alpha.items()}
            if states[0] == "UNSTABLE":
                base_unstable_ids.append(sid)
            for mask, state in states.items():
                portfolio_class_counts[mask][state] += 1
                if mask and state == "UNSTABLE":
                    cardinality_any_unstable[mask.bit_count()] += 1
                if mask and state == "INDETERMINATE":
                    cardinality_any_indeterminate[mask.bit_count()] += 1

            if any(state == "INDETERMINATE" for state in states.values()):
                result["metric_reconciliation"] = "INDETERMINATE_PORTFOLIO_CLASS"
            else:
                blockers = tuple(
                    sorted(
                        (
                            mask
                            for mask in range(1, 16)
                            if states[mask] == "UNSTABLE"
                            and all(
                                states[submask] == "STABLE"
                                for submask in range(16)
                                if submask != mask and (submask & mask) == submask
                            )
                        ),
                        key=lambda mask: (mask.bit_count(), mask),
                    )
                )
                x_value = max(alpha[mask] for mask in range(16) if mask != 15)
                y_value = alpha[15]
                margin = min(y_value, -x_value)
                exact_label = exact_blocker_label(blockers)
                exact_set_counts[exact_label] += 1
                for mask in blockers:
                    blocker_membership[mask] += 1
                    if mask.bit_count() == 3:
                        triple_membership[mask] += 1
                kappa = min((mask.bit_count() for mask in blockers), default=None)
                kappa_counts[kappa] += 1

                metric_masks = mask_list(row["minimal_unstable_blocker_masks"])
                if tuple(sorted(blockers)) != metric_masks:
                    metric_mask_mismatches.append(sid)
                if not (
                    np.isclose(x_value, float(row["x_s_max_proper_alpha"]), rtol=1e-12, atol=1e-14)
                    and np.isclose(y_value, float(row["y_s_h4_alpha_original"]), rtol=1e-12, atol=1e-14)
                ):
                    metric_xy_mismatches.append(sid)

                result.update(
                    {
                        "x_s_max_proper_alpha_s-1": x_value,
                        "y_s_H4_alpha_s-1": y_value,
                        "m_comp_min_y_minus_x_s-1": margin,
                        "minimal_blocker_masks": ";".join(map(str, blockers)) if blockers else "NONE",
                        "minimal_blocker_ids": " | ".join(mask_name(mask) for mask in blockers) if blockers else "NONE",
                        "exact_blocker_set_category": exact_label,
                        "n_minimal_blockers": len(blockers),
                        "kappa_minimal_blocker_order": kappa,
                        "h4_unstable": "TRUE" if states[15] == "UNSTABLE" else "FALSE",
                        "metric_reconciliation": "PASS" if tuple(sorted(blockers)) == metric_masks else "MISMATCH",
                    }
                )

                scenario_results[sid] = {
                    "scenario_id": sid,
                    "scenario_index_zero_based": int(row["scenario_index_zero_based"]),
                    "valid": True,
                    "blockers": blockers,
                    "label": exact_label,
                    "kappa": kappa,
                    "x": x_value,
                    "y": y_value,
                    "margin": margin,
                    "alpha": alpha,
                    "load_scale": float(row["load_scale"]),
                    "epsilon_l2": float(row["epsilon_l2"]),
                    "saturation_count": int(row["saturation_count"]),
                }

        records.append(result)
        provenance.append(
            {
                "scenario_id": sid,
                "scenario_valid": bool(is_valid),
                "source_run_id": SOURCE_RUN_ID,
                "MC_CASES_source_data_line_numbers_by_mask_ascending": source_lines,
                "SCENARIO_METRICS_source_data_line_number": int(row["source_row_line"]),
                "MC_CASES_source_path": MC_CASES.relative_to(REPO).as_posix(),
                "SCENARIO_METRICS_source_path": SCENARIO_METRICS.relative_to(REPO).as_posix(),
                "included_in_dual_boundary_plane": bool(is_valid),
                "included_in_blocker_set_flow": bool(is_valid),
                "figure_ids": "IAS26-111-F01;IAS26-111-F02" if is_valid else "NONE_PHYSICALLY_INFEASIBLE",
            }
        )

    valid_results = [scenario_results[sid] for sid in scenario_results]
    n_valid = len(valid_results)
    exact_categories = {
        "none": exact_set_counts.get("NONE", 0),
        "h4_only": exact_set_counts.get("H4_ONLY", 0),
        "alternative": sum(count for label, count in exact_set_counts.items() if label not in {"NONE", "H4_ONLY"}),
        "valid": n_valid,
    }
    total_h4_unstable = sum(1 for record in valid_results if record["y"] > TAU_DEC)
    alternative_witness_cases = [record for record in valid_results if record["label"] not in {"NONE", "H4_ONLY"}]
    singleton_masks = [mask for mask in range(1, 16) if mask.bit_count() == 1]
    pair_masks = [mask for mask in range(1, 16) if mask.bit_count() == 2]
    singleton_unstable_cells = sum(portfolio_class_counts[m]["UNSTABLE"] for m in singleton_masks)
    pair_unstable_cells = sum(portfolio_class_counts[m]["UNSTABLE"] for m in pair_masks)
    singleton_indeterminate_cells = sum(portfolio_class_counts[m]["INDETERMINATE"] for m in singleton_masks)
    pair_indeterminate_cells = sum(portfolio_class_counts[m]["INDETERMINATE"] for m in pair_masks)
    triple_unstable_counts = {mask: portfolio_class_counts[mask]["UNSTABLE"] for mask in (7, 11, 13, 14)}
    exact_pattern_expected = {
        "NONE": 412,
        "H4_ONLY": 190,
        "{30,33,35}": 77,
        "{30,33,35} | {30,33,37}": 151,
        "{30,33,35} | {30,33,37} | {33,35,37}": 137,
    }
    exact_pattern_actual = {label: int(exact_set_counts.get(label, 0)) for label in exact_pattern_expected}
    extra_patterns = {label: count for label, count in exact_set_counts.items() if label not in exact_pattern_expected}

    comparison_fields = {
        "valid_count": n_valid == EXPECTED["valid"],
        "none_count": exact_categories["none"] == EXPECTED["none"],
        "h4_only_count": exact_categories["h4_only"] == EXPECTED["h4_only"],
        "alternative_witness_count": exact_categories["alternative"] == EXPECTED["alternative"],
        "exact_pattern_counts": exact_pattern_actual == exact_pattern_expected and not extra_patterns,
        "all_1000_ids_retained": len(records) == 1000,
        "frozen_manifest_id_set": len(manifest_raw) == 1000 and set(manifest_raw["scenario_id"]) == set(metrics["scenario_id"]),
        "all_portfolios_present_once": not any("each exact V4 mask" in item for item in structural_errors),
        "physical_feasibility_gate_consistency": not physical_mismatch_ids,
        "no_outcome_dependent_filtering": outcome_filter_true_count == 0,
        "valid_rows_have_complete_finite_spectra_and_pass_status": bool(spectra_complete.all() and alpha_finite.all() and valid_cases_pass.all()),
        "frozen_threshold_port_order_and_original_g": frozen_tau == TAU_DEC and frozen_buses == BUSES and frozen_g == 0.03625,
        "no_base_instability": not base_unstable_ids,
        "no_indeterminate_valid_portfolios": all(
            counter["INDETERMINATE"] == 0 for counter in portfolio_class_counts.values()
        ),
        "no_unstable_singletons": singleton_unstable_cells == 0,
        "no_unstable_pairs": pair_unstable_cells == 0,
        "metric_mask_reconciliation": len(metric_mask_mismatches) == 0,
        "metric_xy_reconciliation": len(metric_xy_mismatches) == 0,
        "no_structural_errors": not structural_errors,
    }
    gate_pass = all(comparison_fields.values())
    gate = {
        "ticket": "IAS26-111",
        "run_id": RUN_ID,
        "source_run_id": SOURCE_RUN_ID,
        "status": "PASS_WITH_LIMITATIONS" if gate_pass else "FAIL",
        "post_hoc_exploratory": True,
        "new_power_flow_DAE_or_eigenvalue_solves": 0,
        "tau_dec_s-1": TAU_DEC,
        "valid_scenario_count": n_valid,
        "invalid_scenario_count": len(records) - n_valid,
        "invalid_ids_retained_and_not_replaced": True,
        "outcome_dependent_filtering_applied": False,
        "outcome_dependent_filtering_true_count_in_frozen_metrics": outcome_filter_true_count,
        "physical_validity_mismatch_count": len(physical_mismatch_ids),
        "exact_categories": exact_categories,
        "exact_blocker_set_pattern_counts": exact_pattern_actual,
        "unexpected_blocker_set_patterns": extra_patterns,
        "kappa_counts_including_none": {
            "kappa_1": int(kappa_counts.get(1, 0)),
            "kappa_2": int(kappa_counts.get(2, 0)),
            "kappa_3": int(kappa_counts.get(3, 0)),
            "kappa_4": int(kappa_counts.get(4, 0)),
            "empty_set": int(kappa_counts.get(None, 0)),
        },
        "singletons": {
            "mask_count": len(singleton_masks),
            "unstable_classifications": int(singleton_unstable_cells),
            "indeterminate_classifications": int(singleton_indeterminate_cells),
            "scenarios_with_any_unstable_singleton": int(sum(1 for sid in valid_ids if any(
                classify(float(value)) == "UNSTABLE"
                for value in cases.loc[(cases.scenario_id == sid) & cases.portfolio_mask.isin(singleton_masks), "alpha_perp"]
            ))),
        },
        "pairs": {
            "mask_count": len(pair_masks),
            "unstable_classifications": int(pair_unstable_cells),
            "indeterminate_classifications": int(pair_indeterminate_cells),
            "scenarios_with_any_unstable_pair": int(sum(1 for sid in valid_ids if any(
                classify(float(value)) == "UNSTABLE"
                for value in cases.loc[(cases.scenario_id == sid) & cases.portfolio_mask.isin(pair_masks), "alpha_perp"]
            ))),
        },
        "triple_portfolios": {
            str(mask): {"portfolio": mask_name(mask), "unstable_scenarios": int(triple_unstable_counts[mask]),
                        "minimal_blocker_memberships": int(triple_membership.get(mask, 0))}
            for mask in (7, 11, 13, 14)
        },
        "H4_unstable_scenarios": int(total_h4_unstable),
        "alternative_witness_with_H4_unstable": int(sum(1 for record in alternative_witness_cases if record["y"] > TAU_DEC)),
        "metric_mask_mismatch_count": len(metric_mask_mismatches),
        "metric_xy_mismatch_count": len(metric_xy_mismatches),
        "structural_errors": structural_errors,
        "gate_checks": comparison_fields,
        "source_sha256": source_hashes(),
    }
    return {
        "gate": gate,
        "records": records,
        "provenance": provenance,
        "scenario_results": valid_results,
        "portfolio_class_counts": portfolio_class_counts,
        "blocker_membership": blocker_membership,
        "exact_set_counts": exact_set_counts,
        "kappa_counts": kappa_counts,
        "triple_membership": triple_membership,
        "metric_mask_mismatches": metric_mask_mismatches,
        "metric_xy_mismatches": metric_xy_mismatches,
        "exact_pattern_counts": exact_pattern_actual,
        "source_case_row_count": len(cases),
        "source_metrics_row_count": len(metrics),
        "source_model_hashes": sorted(set(cases["model_hash"].dropna().astype(str))),
        "source_equilibrium_hash_count": int(cases["equilibrium_hash"].nunique(dropna=True)),
    }


def write_gate(analysis: dict[str, object]) -> None:
    dump_json(RUN_ROOT / "reports/IAS26-111_RECONCILIATION_GATE.json", analysis["gate"])


def write_tables_and_claims(analysis: dict[str, object]) -> None:
    derived = RUN_ROOT / "derived"
    claims = RUN_ROOT / "claims"
    reports = RUN_ROOT / "reports"
    for directory in (derived, claims, reports, RUN_ROOT / "provenance"):
        directory.mkdir(parents=True, exist_ok=True)

    atlas = pd.DataFrame(analysis["records"])
    atlas.to_csv(derived / "SCENARIO_ATLAS.csv", index=False, float_format="%.17g")
    pd.DataFrame(analysis["provenance"]).to_csv(RUN_ROOT / "provenance/SCENARIO_SOURCE_LINES.csv", index=False)
    valid_atlas = atlas[atlas["scenario_valid"]].copy()
    summary_rows = []
    for category, group in valid_atlas.groupby("exact_blocker_set_category", sort=False):
        summary_row: dict[str, object] = {"exact_blocker_set_category": category, "scenario_count": len(group)}
        for column, short in (
            ("x_s_max_proper_alpha_s-1", "x_s"),
            ("y_s_H4_alpha_s-1", "y_s"),
            ("m_comp_min_y_minus_x_s-1", "m_comp"),
        ):
            values = group[column].to_numpy(dtype=float)
            for q, label in ((0.05, "q05"), (0.25, "q25"), (0.5, "median"), (0.75, "q75"), (0.95, "q95")):
                summary_row[f"{short}_{label}_s-1"] = float(np.quantile(values, q))
        summary_rows.append(summary_row)
    pd.DataFrame(summary_rows).to_csv(derived / "BOUNDARY_COORDINATE_SUMMARY.csv", index=False, float_format="%.17g")

    portfolio_rows = []
    for mask in MASKS:
        counts = analysis["portfolio_class_counts"][mask]
        portfolio_rows.append(
            {
                "portfolio_mask": mask,
                "portfolio_id": "BASE" if mask == 0 else "+".join(str(b) for bit, b in enumerate(BUSES) if mask & (1 << bit)),
                "portfolio_set": mask_name(mask),
                "cardinality": mask.bit_count(),
                "n_valid_scenarios": EXPECTED["valid"],
                "stable_scenarios": int(counts["STABLE"]),
                "unstable_scenarios": int(counts["UNSTABLE"]),
                "indeterminate_scenarios": int(counts["INDETERMINATE"]),
                "decision_threshold_tau_dec_s-1": TAU_DEC,
            }
        )
    pd.DataFrame(portfolio_rows).to_csv(derived / "PORTFOLIO_CLASS_COUNTS.csv", index=False)

    mask_rows = []
    for mask in range(1, 16):
        mask_rows.append(
            {
                "portfolio_mask": mask,
                "portfolio_set": mask_name(mask),
                "cardinality": mask.bit_count(),
                "minimal_blocker_scenario_memberships": int(analysis["blocker_membership"].get(mask, 0)),
                "scenario_frequency_over_valid_n": analysis["blocker_membership"].get(mask, 0) / EXPECTED["valid"],
                "interpretation": "membership count; scenarios with multiple minimal blockers contribute once per blocker",
            }
        )
    pd.DataFrame(mask_rows).to_csv(derived / "BLOCKER_MASK_COUNTS.csv", index=False)

    set_rows = []
    for label, count in sorted(analysis["exact_set_counts"].items(), key=lambda item: (-item[1], item[0])):
        set_rows.append(
            {
                "exact_minimal_blocker_set": label,
                "scenario_count": int(count),
                "share_valid": count / EXPECTED["valid"],
                "denominator": EXPECTED["valid"],
                "interpretation": "descriptive frequency in the declared synthetic ensemble; not a real-world probability",
            }
        )
    pd.DataFrame(set_rows).to_csv(derived / "BLOCKER_SET_COUNTS.csv", index=False)

    kappa_rows = []
    for kappa in (None, 1, 2, 3, 4):
        count = int(analysis["kappa_counts"].get(kappa, 0))
        kappa_rows.append(
            {
                "kappa": "NONE" if kappa is None else kappa,
                "scenario_count": count,
                "share_valid": count / EXPECTED["valid"],
                "denominator": EXPECTED["valid"],
                "definition": "minimum cardinality among inclusion-minimal unstable V4 masks; NONE when no such mask exists",
            }
        )
    pd.DataFrame(kappa_rows).to_csv(derived / "BLOCKER_ORDER_COUNTS.csv", index=False)

    triple_rows = []
    for mask in (7, 11, 13, 14):
        portfolio_counts = analysis["portfolio_class_counts"][mask]
        triple_rows.append(
            {
                "portfolio_mask": mask,
                "triple": mask_name(mask),
                "unstable_scenarios": int(portfolio_counts["UNSTABLE"]),
                "minimal_blocker_memberships": int(analysis["triple_membership"].get(mask, 0)),
                "indeterminate_scenarios": int(portfolio_counts["INDETERMINATE"]),
                "denominator_valid": EXPECTED["valid"],
            }
        )
    pd.DataFrame(triple_rows).to_csv(derived / "TRIPLE_AUDIT.csv", index=False)

    singleton_masks = [m for m in range(1, 16) if m.bit_count() == 1]
    pair_masks = [m for m in range(1, 16) if m.bit_count() == 2]
    audit_rows = []
    for cardinality, masks in ((1, singleton_masks), (2, pair_masks)):
        unstable = sum(analysis["portfolio_class_counts"][m]["UNSTABLE"] for m in masks)
        indeterminate = sum(analysis["portfolio_class_counts"][m]["INDETERMINATE"] for m in masks)
        audit_rows.append(
            {
                "cardinality": cardinality,
                "portfolio_count": len(masks),
                "scenario_portfolio_cells": len(masks) * EXPECTED["valid"],
                "unstable_cells": int(unstable),
                "indeterminate_cells": int(indeterminate),
                "any_unstable_portfolio": bool(unstable),
                "gate": "PASS" if unstable == 0 and indeterminate == 0 else "FAIL",
                "threshold_tau_dec_s-1": TAU_DEC,
            }
        )
    pd.DataFrame(audit_rows).to_csv(derived / "SINGLETON_PAIR_AUDIT.csv", index=False)

    counterexamples = []
    for record in analysis["scenario_results"]:
        if record["label"] == "H4_ONLY":
            continue
        if record["label"] == "NONE":
            reason = "NO_MINIMAL_UNSTABLE_PORTFOLIO_IN_V4"
        else:
            reason = "H4_UNSTABLE_BUT_NOT_MINIMAL; PROPER_TRIPLE_WITNESS_PRESENT"
        counterexamples.append(
            {
                "scenario_id": record["scenario_id"],
                "counterexample_to": "fixed_H4_minimal_witness_across_all_valid_scenarios",
                "counterexample_class": reason,
                "exact_minimal_blocker_set": record["label"],
                "minimal_blocker_masks": ";".join(map(str, record["blockers"])) if record["blockers"] else "NONE",
                "kappa": record["kappa"],
                "x_s_max_proper_alpha_s-1": record["x"],
                "y_s_H4_alpha_s-1": record["y"],
                "m_comp_s-1": record["margin"],
                "load_scale": record["load_scale"],
                "epsilon_l2": record["epsilon_l2"],
                "saturation_count": record["saturation_count"],
                "source_run_id": SOURCE_RUN_ID,
            }
        )
    pd.DataFrame(counterexamples).to_csv(derived / "COUNTEREXAMPLES.csv", index=False, float_format="%.17g")

    ledger = [
        {"claim_id": "IAS26-111-C1", "claim": "The valid synthetic V4 cohort partitions into 412 with no minimal blocker, 190 H4-only, and 365 with alternative minimal witnesses.", "status": "SUPPORTED_WITHIN_SCOPE", "evidence": "Exact reclassification of frozen MC_CASES; 967 valid IDs; 0 mask mismatches versus frozen SCENARIO_METRICS.", "limitation": "Post hoc; frozen synthetic ensemble; restricted to V4."},
        {"claim_id": "IAS26-111-C2", "claim": "No singleton or pair portfolio is unstable in any of the 967 valid scenarios under tau_dec=1e-8.", "status": "SUPPORTED_WITHIN_SCOPE", "evidence": "SINGLETON_PAIR_AUDIT.csv and PORTFOLIO_CLASS_COUNTS.csv.", "limitation": "Only the 16 portfolios in V4 and the frozen validity gate."},
        {"claim_id": "IAS26-111-C3", "claim": "The minimum blocker order is 3 in 365 scenarios, 4 in 190, and undefined in 412.", "status": "SUPPORTED_WITHIN_SCOPE", "evidence": "BLOCKER_ORDER_COUNTS.csv; exact inclusion-minimal masks in SCENARIO_ATLAS.csv.", "limitation": "V4 only; not a network-wide or V9 blocker-order statement."},
        {"claim_id": "IAS26-111-C4", "claim": "The alternative minimal blockers observed are {30,33,35}, {30,33,37}, and {33,35,37}; every observed minimal triple includes bus 33.", "status": "SUPPORTED_DESCRIPTIVE", "evidence": "BLOCKER_MASK_COUNTS.csv and TRIPLE_AUDIT.csv.", "limitation": "Association only; no causal anchor/backbone claim."},
        {"claim_id": "IAS26-111-C5", "claim": "m_comp=min(y,-x) is a signed coordinate margin in the (x,y) boundary plane.", "status": "SUPPORTED_DEFINITION", "evidence": "SCENARIO_ATLAS.csv and report definition.", "limitation": "Not a physical probability or universal stability margin."},
        {"claim_id": "IAS26-111-C6", "claim": "The collective closure mechanism remains invariant when witness identity changes.", "status": "UNSUPPORTED_NOT_TESTED", "evidence": "No M_ii or Q_B factorization was reconstructed in IAS26-111.", "limitation": "Requires IAS26-112."},
        {"claim_id": "IAS26-111-C7", "claim": "Bus 33 is a causal anchor/backbone for blocker switching.", "status": "UNSUPPORTED_NOT_TESTED", "evidence": "Incidence is descriptive; no singular-vector, ablation, or null-model test was run.", "limitation": "Requires IAS26-113."},
        {"claim_id": "IAS26-111-C8", "claim": "The V4 pattern generalizes to V9 or IEEE-39-wide blockers.", "status": "UNSUPPORTED_NOT_TESTED", "evidence": "Only the 16 V4 portfolios were reanalyzed.", "limitation": "Requires a separately authorized V9 campaign."},
        {"claim_id": "IAS26-111-C9", "claim": "The scenario frequencies are real-world probabilities.", "status": "UNSUPPORTED_INVALID_INTERPRETATION", "evidence": "The input is a declared synthetic stress ensemble.", "limitation": "Report only descriptive ensemble frequencies."},
    ]
    pd.DataFrame(ledger).to_csv(claims / "CLAIM_LEDGER.csv", index=False)

    safe_claims = """# IAS26-111: claims supported within scope\n\nStatus: PASS_WITH_LIMITATIONS. This is a post-hoc descriptive reanalysis of the frozen synthetic IAS26-060 run. No power-flow, DAE, or eigenvalue solve was executed.\n\n- Among 967 physically valid scenarios in the declared synthetic ensemble, 412 had no inclusion-minimal unstable V4 portfolio, 190 had H4={30,33,35,37} as the sole minimal blocker, and 365 had one or more alternative minimal triple blockers.\n- Across the same 967 scenarios, all 4 singleton portfolios and all 6 pair portfolios were stable; there were no indeterminate portfolio classifications at tau_dec=1e-8.\n- The minimum V4 blocker order was 3 for 365 scenarios, 4 for 190, and undefined for 412. This is not a full-network/V9 result.\n- The observed minimal triples were {30,33,35} (365 scenario memberships), {30,33,37} (288), and {33,35,37} (137). All observed minimal triples include bus 33; this is descriptive and not evidence of causality.\n- x_s=max over proper H4 subsets alpha_perp, y_s=alpha_perp(H4), and m_comp=min(y_s,-x_s) are reported per valid scenario. m_comp is a signed coordinate margin in the (x,y) plane, not a probability or a universal physical stability margin.\n- Recomputed minimal masks and boundary coordinates reconcile to the frozen SCENARIO_METRICS table for all 967 valid IDs with zero mismatches. All 1,000 IDs remain visible; 33 physically invalid IDs were not replaced and were not admitted to the valid-cohort counts.\n\nAll rates are descriptive frequencies under this specific synthetic ensemble, not real-world probabilities.\n"""
    unsupported = """# IAS26-111: failed or unsupported claims\n\n## Refuted for the frozen valid cohort\n\n- A fixed H4 minimal witness is not present in every valid scenario: 777/967 valid scenarios do not have H4 as their sole minimal blocker. Of these, 412 have no minimal unstable V4 portfolio and 365 have alternative minimal triple witnesses. This refutes universal persistence of this exact witness within the tested V4 cohort; it does not refute collective instability more generally.\n\n## Not tested by IAS26-111\n\n- Mechanism invariance: no physical local factors I+M_ii or collective I+Q_B matrices were reconstructed.\n- Bus 33 as a causal anchor/backbone: no singular-vector participation, leave-one-port-out ablation, or combinatorial null test was run.\n- IEEE-39-wide or V9 blocker generalization: only V4's 16 portfolios were included.\n- Why the witness family varies with operating point: this ticket is descriptive and does not isolate causes.\n- Any real-world failure probability: the ensemble is synthetic and its frequencies are not naturalistic probabilities.\n- A universal effect of retuning g or an eigenlocus explanation: no new g sweep or Q_H analysis was performed here.\n- Poster-level replacement claims or figures: the poster was not modified or rebuilt.\n"""
    (claims / "SAFE_CLAIMS.md").write_text(safe_claims, encoding="utf-8")
    (claims / "FAILED_OR_UNSUPPORTED_CLAIMS.md").write_text(unsupported, encoding="utf-8")

    gate = analysis["gate"]
    report = f"""# IAS26-111 - Exact V4 blocker atlas\n\n- Run: `{RUN_ID}`\n- Source: `{SOURCE_RUN_ID}`\n- Status: **{gate['status']}**\n- Type: post-hoc descriptive reanalysis; no new power-flow, DAE, or eigenvalue solves.\n- Frozen threshold: `tau_dec = {TAU_DEC:.0e} s^-1`.\n- Valid scenarios: {EXPECTED['valid']}; physically invalid IDs retained without replacement: 33.\n- Outcome-dependent filtering: none.\n\n## Reconciliation gate\n\nThe exact required partition is 412 no-blocker + 190 H4-only + 365 alternative-witness = 967 valid scenarios. The recomputation produced {gate['exact_categories']['none']} + {gate['exact_categories']['h4_only']} + {gate['exact_categories']['alternative']} = {gate['exact_categories']['valid']}. The scenario-level blocker-mask mismatch count against frozen `SCENARIO_METRICS.csv` is {gate['metric_mask_mismatch_count']}; the `(x_s,y_s)` mismatch count is {gate['metric_xy_mismatch_count']}.\n\n## Results\n\n| Exact minimal blocker set | Scenarios | Share of valid cohort |\n|---|---:|---:|\n"""
    for label, count in sorted(analysis["exact_set_counts"].items(), key=lambda item: (-item[1], item[0])):
        report += f"| `{label}` | {count} | {count / EXPECTED['valid']:.2%} |\n"
    report += "\nThe minimum blocker order counts are reported in `../derived/BLOCKER_ORDER_COUNTS.csv`; the complete row-level result is in `../derived/SCENARIO_ATLAS.csv`. For each valid scenario, `x_s=max_{R proper subset H4} alpha_perp(R)`, `y_s=alpha_perp(H4)`, and `m_comp=min(y_s,-x_s)`. The latter is a signed coordinate margin in the `(x,y)` plane, not a probability.\n\n"
    report += f"H4 was unstable in {gate['H4_unstable_scenarios']} valid scenarios. It was minimal in 190 and nonminimal in 365 alternative-witness scenarios. Every alternative witness is represented in `COUNTEREXAMPLES.csv`.\n\n"
    report += "## Scope boundary\n\nThis atlas establishes descriptive witness switching only within the frozen synthetic ensemble and the 16 V4 portfolios. It does not establish mechanism invariance, bus-33 causality, V9/IEEE-39-wide generality, or naturalistic failure probabilities. No poster source, poster figure, or historical run was modified. See `../claims/SAFE_CLAIMS.md` and `../claims/FAILED_OR_UNSUPPORTED_CLAIMS.md`.\n"
    (reports / "IAS26-111_REPORT.md").write_text(report, encoding="utf-8")


def category_label(record: dict[str, object]) -> str:
    return str(record["label"])


def write_figures(analysis: dict[str, object]) -> None:
    figures = RUN_ROOT / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    results = analysis["scenario_results"]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 14,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    labels_order = [
        "NONE",
        "H4_ONLY",
        "{30,33,35}",
        "{30,33,35} | {30,33,37}",
        "{30,33,35} | {30,33,37} | {33,35,37}",
    ]
    labels_full = {
        "NONE": "No minimal blocker (n=412)",
        "H4_ONLY": "H4 only (n=190)",
        "{30,33,35}": "{30,33,35} only (n=77)",
        "{30,33,35} | {30,33,37}": "Two triples (n=151)",
        "{30,33,35} | {30,33,37} | {33,35,37}": "Three triples (n=137)",
    }

    # Figure 1: x-y boundary plane, with the preregistered deadband noted.
    fig, ax = plt.subplots(figsize=(9.5, 7.0), constrained_layout=True)
    xs = np.asarray([float(record["x"]) for record in results])
    ys = np.asarray([float(record["y"]) for record in results])
    xpad = max(0.025, 0.07 * (xs.max() - xs.min()))
    ypad = max(0.025, 0.07 * (ys.max() - ys.min()))
    xmin, xmax = xs.min() - xpad, xs.max() + xpad
    ymin, ymax = ys.min() - ypad, ys.max() + ypad
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.add_patch(
        Rectangle(
            (xmin, 0),
            -xmin,
            ymax,
            facecolor="#D8EDE6",
            edgecolor="none",
            alpha=0.44,
            zorder=0,
        )
    )
    ax.axvline(0, color="#303942", linewidth=1.05, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="#303942", linewidth=1.05, linestyle=(0, (4, 3)), zorder=1)
    for label in labels_order:
        group = [record for record in results if category_label(record) == label]
        if not group:
            continue
        ax.scatter(
            [record["x"] for record in group],
            [record["y"] for record in group],
            s=25,
            color=COLORS[label],
            edgecolors="white",
            linewidths=0.3,
            alpha=0.78,
            label=labels_full[label],
            rasterized=False,
            zorder=2,
        )
    ax.text(xmin + 0.02 * (xmax - xmin), ymax - 0.04 * (ymax - ymin), "H4-minimal quadrant", color="#376B5E", va="top", fontsize=9, weight="bold")
    ax.set_xlabel(r"$x_s = \max_{R\subsetneq H4}\,\alpha_{\perp,s}(R)$ (s$^{-1}$)")
    ax.set_ylabel(r"$y_s = \alpha_{\perp,s}(H4)$ (s$^{-1}$)")
    ax.set_title("V4 stability-boundary coordinates by exact minimal-blocker set", loc="left", weight="bold", pad=12)
    ax.text(
        0.0,
        -0.16,
        "Dashed axes show x=0 and y=0. Decisions use the frozen ±1e-8 s^-1 deadband (too narrow to resolve at this scale). n=967 valid IDs.",
        transform=ax.transAxes,
        fontsize=8,
        color="#505A64",
    )
    ax.legend(loc="lower right", frameon=True, framealpha=0.96, title="Exact minimal-blocker set", title_fontsize=8)
    fig.savefig(figures / "IAS26-111_F01_DUAL_BOUNDARY_PLANE.pdf", bbox_inches="tight", metadata={"Title": "IAS26-111 V4 dual-boundary plane", "Subject": "Post-hoc descriptive frozen-scenario reanalysis"})
    fig.savefig(figures / "IAS26-111_F01_DUAL_BOUNDARY_PLANE.svg", bbox_inches="tight", metadata={"Title": "IAS26-111 V4 dual-boundary plane"})
    plt.close(fig)

    # Figure 2: categorical partition rendered as a flow; not a trajectory.
    counts = analysis["exact_pattern_counts"]
    flow_labels = [
        ("NONE", "No minimal blocker in V4"),
        ("H4_ONLY", "H4={30,33,35,37} only"),
        ("{30,33,35}", "{30,33,35} only"),
        ("{30,33,35} | {30,33,37}", "{30,33,35} + {30,33,37}"),
        ("{30,33,35} | {30,33,37} | {33,35,37}", "Three minimal triples"),
    ]
    fig, ax = plt.subplots(figsize=(12.4, 7.8), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.98, "Observed V4 minimal-blocker-set partition", fontsize=15, weight="bold", va="top", color="#202932")
    ax.text(0.02, 0.925, "967 valid scenarios; categorical partition of IDs, not a temporal or causal flow", fontsize=9.5, va="top", color="#58636E")

    total_top, total_bottom = 0.86, 0.12
    total_height = total_top - total_bottom
    left_x0, left_x1 = 0.06, 0.19
    mid_x, right_x0, right_x1 = 0.52, 0.74, 0.98
    ax.add_patch(Rectangle((left_x0, total_bottom), left_x1 - left_x0, total_height, facecolor="#E4E9ED", edgecolor="#4F5D69", linewidth=1.0))
    ax.text((left_x0 + left_x1) / 2, 0.5, "Valid V4\nscenarios\nn=967", ha="center", va="center", fontsize=10, weight="bold", color="#27333D")

    y_cursor = total_top
    for key, display in flow_labels:
        count = int(counts.get(key, 0))
        height = total_height * count / EXPECTED["valid"]
        y0 = y_cursor - height
        y1 = y_cursor
        color = COLORS[key]
        # Smooth ribbon, retaining the count-proportional band thickness.
        verts = [
            (left_x1, y0), (left_x1 + 0.15, y0), (mid_x - 0.10, y0), (mid_x, y0),
            (mid_x, y1), (mid_x - 0.10, y1), (left_x1 + 0.15, y1), (left_x1, y1),
            (left_x1, y0),
        ]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.LINETO,
                 MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="white", linewidth=0.7, alpha=0.72, zorder=1))
        ax.add_patch(Rectangle((right_x0, y0), right_x1 - right_x0, height, facecolor=color, edgecolor="white", linewidth=0.8, zorder=2))
        percentage = count / EXPECTED["valid"]
        ax.text(right_x0 + 0.012, (y0 + y1) / 2, f"{display}\n{count}  ({percentage:.2%})", ha="left", va="center", fontsize=8.2, color="white" if key not in {"NONE", "H4_ONLY"} else "#17212A", weight="bold", zorder=3)
        y_cursor = y0

    ax.text(left_x0, 0.075, "One row per frozen scenario ID", fontsize=8, color="#596570")
    ax.text(right_x0, 0.075, "Counts sum to 967; invalid IDs are outside this partition", fontsize=8, color="#596570")
    fig.savefig(figures / "IAS26-111_F02_BLOCKER_SET_FLOW.pdf", bbox_inches="tight", metadata={"Title": "IAS26-111 blocker-set partition flow", "Subject": "Categorical partition, not temporal or causal flow"})
    fig.savefig(figures / "IAS26-111_F02_BLOCKER_SET_FLOW.svg", bbox_inches="tight", metadata={"Title": "IAS26-111 blocker-set partition flow"})
    plt.close(fig)


def finalize_manifest(analysis: dict[str, object]) -> None:
    if analysis["gate"]["status"] != "PASS_WITH_LIMITATIONS":
        raise RuntimeError("Cannot finalize a failed reconciliation gate")
    current_sources = source_hashes()
    if current_sources != analysis["gate"]["source_sha256"]:
        raise RuntimeError("A frozen input hash changed during IAS26-111")
    output_files = sorted(
        path for path in RUN_ROOT.rglob("*")
        if path.is_file() and path.name not in {"manifest.json", "SHA256SUMS.txt"}
    )
    output_hash_map = {path.relative_to(RUN_ROOT).as_posix(): sha256(path) for path in output_files}
    manifest = {
        "ticket": "IAS26-111",
        "run_id": RUN_ID,
        "status": "PASS_WITH_LIMITATIONS",
        "created_utc": RUN_ID[:16],
        "post_hoc_exploratory": True,
        "source_run_id": SOURCE_RUN_ID,
        "source_paths": {
            "MC_CASES": MC_CASES.relative_to(REPO).as_posix(),
            "SCENARIO_METRICS": SCENARIO_METRICS.relative_to(REPO).as_posix(),
            "frozen_config": CONFIG.relative_to(REPO).as_posix(),
            "frozen_scenario_manifest": SCENARIOS.relative_to(REPO).as_posix(),
            "prior_read_only_atlas_memo": SOURCE_MEMO.relative_to(REPO).as_posix(),
        },
        "source_sha256": current_sources,
        "source_census": {
            "MC_CASES_rows_total": 17000,
            "original_treatment_cases_used": analysis["source_case_row_count"],
            "SCENARIO_METRICS_rows": analysis["source_metrics_row_count"],
            "model_hashes": analysis["source_model_hashes"],
            "distinct_equilibrium_hashes_across_cases": analysis["source_equilibrium_hash_count"],
        },
        "analysis_contract": {
            "V4_bus_order": list(BUSES),
            "mask_bit_order": "bit 0..3 maps to bus 30,33,35,37",
            "tau_dec_s-1": TAU_DEC,
            "stable": "alpha_perp < -tau_dec",
            "unstable": "alpha_perp > tau_dec",
            "indeterminate": "-tau_dec <= alpha_perp <= tau_dec",
            "minimal_unstable": "unstable nonempty mask for which every strict subset including BASE is stable",
            "x_s": "maximum alpha_perp across the 15 strict proper subsets of H4, including BASE",
            "y_s": "alpha_perp of mask 15 / H4",
            "m_comp": "min(y_s,-x_s); signed coordinate margin in (x,y) plane, not a probability",
            "validity": "scenario_valid from frozen SCENARIO_METRICS; invalid IDs retained, not replaced, and excluded from atlas frequencies",
            "new_solves": 0,
            "poster_modified": False,
        },
        "gate": analysis["gate"],
        "output_sha256": output_hash_map,
        "environment": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "scope_limitations": [
            "Post-hoc descriptive analysis of one frozen synthetic scenario ensemble.",
            "V4 only; no V9/general IEEE-39 blocker census.",
            "No Q factorization, local-factor audit, causal anchor test, g sweep, or mode-identity resolution.",
            "Scenario frequencies are not real-world probabilities.",
        ],
    }
    dump_json(RUN_ROOT / "manifest.json", manifest)
    checksum_files = sorted(path for path in RUN_ROOT.rglob("*") if path.is_file() and path.name != "SHA256SUMS.txt")
    checksum_text = "".join(f"{sha256(path)}  {path.relative_to(RUN_ROOT).as_posix()}\n" for path in checksum_files)
    (RUN_ROOT / "SHA256SUMS.txt").write_text(checksum_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("audit", "tables", "figures", "finalize"), required=True)
    args = parser.parse_args()
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    analysis = load_and_analyze()
    if args.phase == "audit":
        write_gate(analysis)
        print(json.dumps(analysis["gate"], indent=2))
        return 0 if analysis["gate"]["status"] == "PASS_WITH_LIMITATIONS" else 2
    if analysis["gate"]["status"] != "PASS_WITH_LIMITATIONS":
        write_gate(analysis)
        print(json.dumps(analysis["gate"], indent=2))
        return 2
    if args.phase == "tables":
        write_tables_and_claims(analysis)
    elif args.phase == "figures":
        write_figures(analysis)
    elif args.phase == "finalize":
        finalize_manifest(analysis)
    print(json.dumps({"run_id": RUN_ID, "phase": args.phase, "gate": analysis["gate"]["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
