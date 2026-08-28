from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PREREG_NAMES = (
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.json", "E05D_HOLDOUT_FREEZE.json",
    "E05D_GRID_STRESS_BRANCH_FREEZE.json", "E05D_STRESS_RULE_FREEZE.json",
    "E05D_CRITICAL_MODE_RULE_FREEZE.json", "E05D_TRACKING_FREEZE.json",
    "E05D_C3C_GATE_FREEZE.json", "E05D_NUMERICAL_FREEZE.json",
)
REPORT_NAMES = (
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.md", "E05D_HOLDOUT_FREEZE.md",
    "E05D_GRID_STRESS_BRANCH_FREEZE.md", "E05D_GRID_STRENGTH_IMPLEMENTATION.md",
    "E05D_STRESS_ENDPOINTS.md", "E05D_CRITICAL_MODE_RULE.md", "E05D_MODE_TRACKING_AUDIT.md",
    "E05D_CONNECTED_POLE_RESULTS.md", "E05D_MARGIN_CONSUMPTION.md", "E05D_COMPOSITION_FAILURES.md",
    "E05D_CONDITIONING_ANALYSIS.md", "E05D_FAILURE_LEDGER.md", "E05D_C3C_DECISION.md",
    "E05D_REVIEWER1_AUDIT.md", "E05D_REVIEWER2_AUDIT.md", "E05D_FINAL_REVIEW.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    artifact = root / "artifacts" / "tx3" / "E05D_weak_grid"
    prereg = artifact / "preregistration"

    def require(condition: bool, message: str) -> None:
        if not condition: errors.append(message)

    hashes_path = prereg / "E05D_PREREGISTRATION_SHA256.json"
    require(hashes_path.is_file(), "missing preregistration hash manifest")
    if not hashes_path.is_file(): return errors
    hashes = json.loads(hashes_path.read_text(encoding="utf-8"))["sha256"]
    for name in PREREG_NAMES:
        path = prereg / name
        require(path.is_file(), f"missing preregistration/{name}")
        if path.is_file(): require(sha256(path) == hashes.get(name), f"preregistration changed: {name}")
    numerical = json.loads((prereg / "E05D_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    for relative, expected in numerical["source_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen source: {relative}")
        if path.is_file(): require(sha256(path) == expected, f"frozen source changed: {relative}")
    candidates = json.loads((prereg / "E05D_CANDIDATE_FREEZE.json").read_text(encoding="utf-8"))
    branches = json.loads((prereg / "E05D_GRID_STRESS_BRANCH_FREEZE.json").read_text(encoding="utf-8"))
    holdout = json.loads((prereg / "E05D_HOLDOUT_FREEZE.json").read_text(encoding="utf-8"))
    require(candidates["candidate_count"] == 15 and candidates["pair_count"] == 10 and candidates["triple_count"] == 5, "candidate freeze changed")
    require([record["line_id"] for record in branches["branches"] if record["included"]] == [f"Line_{index}" for index in range(1, 35)], "branch freeze changed")
    require(holdout["seed"] == 20260830 and len(holdout["operating_points"]) == 8, "holdout freeze changed")
    tables = artifact / "tables"
    holdout_dir = artifact / "holdout"
    limits = pd.read_parquet(tables / "e05d_seed_stress_limits.parquet")
    grid = pd.read_parquet(tables / "e05d_grid_strength.parquet")
    ledger = pd.read_parquet(holdout_dir / "e05d_failure_ledger.parquet")
    require(len(limits) == 8 and limits["status"].eq("SUCCESS").all(), "8/8 physical/eigen baselines not retained")
    require(limits["calibration_status"].eq("BASELINE_INELIGIBLE_DAMPING_BELOW_0.05").all(), "eligibility classification changed")
    require((limits["minimum_oscillatory_damping_ratio"] < 0.05).all(), "a baseline is not below the frozen screen")
    require(limits["kappa_end"].isna().all(), "an inadmissible kappa_end was manufactured")
    require(len(grid) == 8 and grid["tau"].eq(0.0).all() and grid["kappa_grid"].eq(1.0).all(), "grid-strength baseline table changed")
    require(len(ledger) == 8 and ledger["status"].eq("SUCCESS").all(), "baseline ledger changed")
    for name in (
        "e05d_critical_modes.parquet", "e05d_mode_tracking.parquet", "e05d_subset_poles.parquet",
        "e05d_connected_externalities.parquet", "e05d_margin_consumption.parquet",
        "e05d_composition_failures.parquet", "e05d_conditioning.parquet",
    ):
        path = tables / name if name == "e05d_critical_modes.parquet" else holdout_dir / name
        require(path.is_file(), f"missing explicit empty table {name}")
        if path.is_file(): require(len(pd.read_parquet(path)) == 0, f"coalition result exists despite baseline stop: {name}")
    final = json.loads((artifact / "E05D_FINAL_DECISION.json").read_text(encoding="utf-8"))
    require(final["terminal_status"] == "BASELINE_INELIGIBLE", "terminal status changed")
    require(final["C3c_WEAK_GRID_MATERIAL_EXTERNALITY"] == "UNRESOLVED", "C3c must remain unresolved")
    require(final["coalition_results_computed"] is False and final["coalition_result_rows"] == 0, "coalition result was computed")
    require(final["materiality_search_closed"] is True and final["c3_mechanism_localization_authorized"] is False, "closure/authorization mismatch")
    require(final["e06_consideration_authorized"] is False and final["E06_executed"] is False and final["E05E_authorized"] is False, "forbidden follow-on stage authorized")
    for name in REPORT_NAMES: require((artifact / name).is_file(), f"missing report {name}")
    stems = (
        "weak_grid_stress", "empty_critical_margin", "all_candidate_rho_heatmap", "lower_order_vs_full_pole",
        "damping_screen", "grid_strength_vs_connected_relevance", "pair_vs_triple", "conditioning",
    )
    for index, stem in enumerate(stems, start=1):
        for suffix in (".pdf", ".png"):
            path = artifact / "figures" / f"figure_E05D_{index}_{stem}{suffix}"
            require(path.is_file() and path.stat().st_size > 1000, f"missing/empty figure {path.name}")
    claims = root / "reports/papers/tx3_connected_intervention_calculus/CLAIMS_EVIDENCE_MATRIX_E05D.md"
    require(claims.is_file() and "| C3c |" in claims.read_text(encoding="utf-8"), "E05D claim matrix missing")
    require(not (root / "artifacts/tx3/E06").exists(), "E06 artifact directory exists")
    manifest = artifact / "manifests" / "PACKAGE_MANIFEST.json"
    if manifest.exists():
        for record in json.loads(manifest.read_text(encoding="utf-8"))["files"]:
            path = root / record["path"]
            require(path.is_file(), f"package file missing: {record['path']}")
            if path.is_file():
                require(path.stat().st_size == record["bytes"], f"package size mismatch: {record['path']}")
                require(sha256(path) == record["sha256"], f"package hash mismatch: {record['path']}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(); errors = validate(args.root.resolve())
    print(json.dumps({"root": str(args.root.resolve()), "status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
