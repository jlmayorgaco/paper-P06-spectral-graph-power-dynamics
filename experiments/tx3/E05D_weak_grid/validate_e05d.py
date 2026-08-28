from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
sys.path.insert(0, str(E05C))
from margin_conditioned import parse_subset_key, union_vertices  # noqa: E402


PREREG_FILES = {
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.json", "E05D_HOLDOUT_FREEZE.json",
    "E05D_GRID_STRESS_BRANCH_FREEZE.json", "E05D_STRESS_RULE_FREEZE.json",
    "E05D_CRITICAL_MODE_RULE_FREEZE.json", "E05D_TRACKING_FREEZE.json",
    "E05D_C3C_GATE_FREEZE.json", "E05D_NUMERICAL_FREEZE.json",
}
REQUIRED_REPORTS = {
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.md", "E05D_HOLDOUT_FREEZE.md",
    "E05D_GRID_STRESS_BRANCH_FREEZE.md", "E05D_GRID_STRENGTH_IMPLEMENTATION.md",
    "E05D_STRESS_ENDPOINTS.md", "E05D_CRITICAL_MODE_RULE.md", "E05D_MODE_TRACKING_AUDIT.md",
    "E05D_CONNECTED_POLE_RESULTS.md", "E05D_MARGIN_CONSUMPTION.md", "E05D_COMPOSITION_FAILURES.md",
    "E05D_CONDITIONING_ANALYSIS.md", "E05D_FAILURE_LEDGER.md", "E05D_C3C_DECISION.md",
    "E05D_REVIEWER1_AUDIT.md", "E05D_REVIEWER2_AUDIT.md", "E05D_FINAL_REVIEW.md",
}
VALID_STATUSES = {
    "SUCCESS", "PF_FAIL", "INIT_FAIL", "EIG_FAIL", "MODE_TRACK_FAIL", "EIGENPAIR_RESIDUAL_FAIL",
    "NONSMOOTH_REGIME", "LIMITER_ACTIVE", "VOLTAGE_LIMIT", "GENERATOR_LIMIT", "BASELINE_UNSTABLE",
    "LOWER_ORDER_UNSTABLE", "NUMERICAL_UNRESOLVED",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    artifact = root / "artifacts" / "tx3" / "E05D_weak_grid"
    prereg = artifact / "preregistration"

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    for name in PREREG_FILES:
        require((prereg / name).is_file(), f"missing preregistration/{name}")
    for name in (
        "E05D_PREREGISTRATION_SHA256.json", "E05D_STRESS_ENDPOINT_FREEZE.json",
        "E05D_CRITICAL_MODE_FREEZE.json", "E05D_CRITICAL_MODE_EIGENVECTORS.npz",
        "e05d_candidates.parquet", "e05d_holdout_manifest.parquet", "e05d_empty_mode_tracking_freeze.parquet",
    ):
        require((prereg / name).is_file(), f"missing preregistration/{name}")
    if errors:
        return errors
    hashes = json.loads((prereg / "E05D_PREREGISTRATION_SHA256.json").read_text(encoding="utf-8"))
    for name in PREREG_FILES:
        require(hashes["sha256"].get(name) == sha256(prereg / name), f"preregistration hash mismatch: {name}")
    numerical = json.loads((prereg / "E05D_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    for relative, expected in numerical["source_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen source {relative}")
        if path.is_file(): require(sha256(path) == expected, f"frozen source hash mismatch: {relative}")
    candidates_freeze = json.loads((prereg / "E05D_CANDIDATE_FREEZE.json").read_text(encoding="utf-8"))
    holdout = json.loads((prereg / "E05D_HOLDOUT_FREEZE.json").read_text(encoding="utf-8"))
    branches = json.loads((prereg / "E05D_GRID_STRESS_BRANCH_FREEZE.json").read_text(encoding="utf-8"))
    endpoint = json.loads((prereg / "E05D_STRESS_ENDPOINT_FREEZE.json").read_text(encoding="utf-8"))
    mode_freeze = json.loads((prereg / "E05D_CRITICAL_MODE_FREEZE.json").read_text(encoding="utf-8"))
    require(candidates_freeze["candidate_count"] == 15 and candidates_freeze["pair_count"] == 10 and candidates_freeze["triple_count"] == 5, "candidate population changed")
    require(holdout["seed"] == 20260830 and len(holdout["operating_points"]) == 8 and not holdout["old_points_reused"], "holdout freeze changed")
    included = [record["line_id"] for record in branches["branches"] if record["included"]]
    require(included == [f"Line_{index}" for index in range(1, 35)], "weak-grid branch set changed")
    require(len(branches["branches"]) == 46 and branches["included_count"] == 34 and branches["excluded_count"] == 12, "branch census changed")
    require(endpoint["coalition_outcomes_inspected_before_freeze"] is False, "endpoint freeze admits coalition inspection")
    for key in ("limits", "trace"):
        path = root / endpoint[f"{key}_path"]
        require(sha256(path) == endpoint[f"{key}_sha256"], f"endpoint {key} artifact changed")
    require(mode_freeze["coalition_outcomes_inspected_before_freeze"] is False, "critical-mode freeze admits coalition inspection")
    for key in ("critical_mode_table", "empty_tracking_table", "eigenvector_archive"):
        path = root / mode_freeze[key]
        require(sha256(path) == mode_freeze[f"{key}_sha256"], f"critical-mode {key} changed")

    tables = artifact / "tables"
    holdout_dir = artifact / "holdout"
    required_tables = (
        tables / "e05d_seed_stress_limits.parquet", tables / "e05d_grid_strength.parquet",
        tables / "e05d_critical_modes.parquet", holdout_dir / "e05d_mode_tracking.parquet",
        holdout_dir / "e05d_subset_poles.parquet", holdout_dir / "e05d_connected_externalities.parquet",
        holdout_dir / "e05d_margin_consumption.parquet", holdout_dir / "e05d_composition_failures.parquet",
        holdout_dir / "e05d_conditioning.parquet", holdout_dir / "e05d_failure_ledger.parquet",
    )
    for path in required_tables: require(path.is_file(), f"missing table {path.name}")
    if any(not path.is_file() for path in required_tables): return errors
    limits = pd.read_parquet(tables / "e05d_seed_stress_limits.parquet")
    grid = pd.read_parquet(tables / "e05d_grid_strength.parquet")
    critical = pd.read_parquet(tables / "e05d_critical_modes.parquet")
    tracks = pd.read_parquet(holdout_dir / "e05d_mode_tracking.parquet")
    results = pd.read_parquet(holdout_dir / "e05d_connected_externalities.parquet")
    builds = pd.read_parquet(holdout_dir / "e05d_failure_ledger.parquet")
    decisions = pd.read_parquet(tables / "e05d_candidate_decisions.parquet")
    final = json.loads((artifact / "E05D_FINAL_DECISION.json").read_text(encoding="utf-8"))
    candidates = candidates_freeze["candidates"]
    union_count = len(union_vertices([parse_subset_key(record["coalition_key"]) for record in candidates]))
    require(len(limits) == 8 and limits["status"].eq("SUCCESS").all(), "baseline calibration incomplete")
    require(len(grid) == 72 and grid["status"].eq("SUCCESS").all(), "grid-strength table incomplete")
    require((grid.groupby("seed_id")["min_SCR"].first() > grid.groupby("seed_id")["min_SCR"].last()).all(), "SCR does not weaken monotonically")
    require(len(critical) == 8 and critical["status"].eq("SUCCESS").all(), "critical-mode freeze incomplete")
    require(len(tracks) == 8 * 9 * union_count, "mode-tracking population incomplete")
    successful_tracks = tracks[tracks["status"] == "SUCCESS"]
    require(float(successful_tracks["BMAC"].min()) >= 0.90, "successful mode track violates BMAC freeze")
    require(float(successful_tracks["eigenpair_residual"].max()) <= 1e-10, "successful mode track violates residual freeze")
    require(len(results) == 1080 and results["seed_id"].nunique() == 8 and results["tau"].nunique() == 9, "connected result population incomplete")
    require(set(results["status"].dropna()).issubset(VALID_STATUSES), "unregistered connected status")
    require(set(builds["status"].dropna()).issubset(VALID_STATUSES), "unregistered build status")
    valid = results[results["status"] == "SUCCESS"]
    require(max(float(valid["composition_closure_real"].abs().max()), float(valid["composition_closure_imag"].abs().max())) <= 1e-10, "complex composition closure failed")
    require(len(decisions) == 15, "candidate decision table incomplete")
    expected_c3c = "SUPPORTED" if decisions["C3c"].eq("SUPPORTED").any() else "REJECTED"
    require(final["C3c_WEAK_GRID_MATERIAL_EXTERNALITY"] == expected_c3c, "C3c aggregation mismatch")
    require(final["c3_mechanism_localization_authorized"] == (expected_c3c == "SUPPORTED"), "mechanism authorization mismatch")
    require(final["e06_consideration_authorized"] is False and final["E06_executed"] is False, "E06 was authorized or executed")
    require(final["E05E_authorized"] is False and final["new_paraemt_scientific_runs"] == 0, "forbidden follow-on stage recorded")
    for name in REQUIRED_REPORTS: require((artifact / name).is_file(), f"missing report {name}")
    for index, stem in enumerate((
        "weak_grid_stress", "empty_critical_margin", "all_candidate_rho_heatmap", "lower_order_vs_full_pole",
        "damping_screen", "grid_strength_vs_connected_relevance", "pair_vs_triple", "conditioning",
    ), start=1):
        for suffix in (".pdf", ".png"):
            path = artifact / "figures" / f"figure_E05D_{index}_{stem}{suffix}"
            require(path.is_file() and path.stat().st_size > 1000, f"missing/empty figure {path.name}")
    claims = root / "reports" / "papers" / "tx3_connected_intervention_calculus" / "CLAIMS_EVIDENCE_MATRIX_E05D.md"
    require(claims.is_file(), "E05D versioned claim matrix missing")
    if claims.is_file():
        text = claims.read_text(encoding="utf-8")
        for marker in ("| C1 |", "| C2 |", "| C3 |", "| C3a |", "| C3b-A |", "| C3b-B |", "| C3c |"):
            require(marker in text, f"claim matrix lacks {marker}")
    e05c = json.loads((root / "artifacts/tx3/E05C_margin_conditioned/E05C_FINAL_DECISION.json").read_text(encoding="utf-8"))
    require(e05c["C3b_A_stress_conditioned_amplification"] == "REJECTED" and e05c["C3b_B_material_margin_consequence"] == "REJECTED", "historical E05C decision changed")
    require(not (root / "artifacts/tx3/E06").exists(), "E06 artifact directory exists")
    package_manifest = artifact / "manifests" / "PACKAGE_MANIFEST.json"
    if package_manifest.exists():
        payload = json.loads(package_manifest.read_text(encoding="utf-8"))
        for record in payload["files"]:
            path = root / record["path"]
            require(path.is_file(), f"package file missing: {record['path']}")
            if path.is_file():
                require(path.stat().st_size == record["bytes"], f"package size mismatch: {record['path']}")
                require(sha256(path) == record["sha256"], f"package hash mismatch: {record['path']}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    errors = validate(args.root.resolve())
    print(json.dumps({"root": str(args.root.resolve()), "status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
