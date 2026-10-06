from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
REQUIRED_REPORTS = {
    "E05C_PREREGISTRATION.md", "E05C_CANDIDATE_FREEZE.md", "E05C_MODE_FAMILY_FREEZE.md",
    "E05C_STRESS_PATH.md", "E05C_MODE_TRACKING_AUDIT.md", "E05C_CONNECTED_POLE_EXTERNALITIES.md",
    "E05C_MARGIN_CONSUMPTION.md", "E05C_COMPOSITION_FAILURES.md", "E05C_CONDITIONING_ANALYSIS.md",
    "E05C_C3B_DECISION.md", "E05C_FAILURE_LEDGER.md", "E05C_REVIEWER1_AUDIT.md",
    "E05C_REVIEWER2_AUDIT.md", "E05C_FINAL_REVIEW.md", "E05C_CLAIMS_EVIDENCE_ADDENDUM.md",
    "INDEPENDENT_REVIEW_README.md",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    artifact = root / "artifacts" / "tx3" / "E05C_margin_conditioned"
    prereg = artifact / "preregistration"

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    for name in (
        "E05C_BASE_HOLDOUT_FREEZE.json", "E05C_BASE_HOLDOUT_MANIFEST.parquet",
        "E05C_CANDIDATE_FREEZE.json", "E05C_MODE_FAMILY_FREEZE.json", "E05C_STRESS_RULE_FREEZE.json",
        "E05C_STRESS_ENDPOINT_FREEZE.json", "E05C_NUMERICAL_FREEZE.json", "E05C_GATE_FREEZE.json",
        "E05C_FIGURE_PANEL_FREEZE.json", "E05C_PREREGISTRATION.md", "e05c_candidates.parquet",
        "e05c_mode_family_freeze.parquet",
    ):
        require((prereg / name).is_file(), f"missing preregistration/{name}")
    if errors:
        return errors
    numerical = json.loads((prereg / "E05C_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    candidate_freeze = json.loads((prereg / "E05C_CANDIDATE_FREEZE.json").read_text(encoding="utf-8"))
    mode_freeze = json.loads((prereg / "E05C_MODE_FAMILY_FREEZE.json").read_text(encoding="utf-8"))
    endpoint_freeze = json.loads((prereg / "E05C_STRESS_ENDPOINT_FREEZE.json").read_text(encoding="utf-8"))
    final = json.loads((artifact / "E05C_FINAL_DECISION.json").read_text(encoding="utf-8"))
    for relative, expected in numerical["source_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen source {relative}")
        if path.is_file():
            require(sha256(path) == expected, f"frozen source hash mismatch: {relative}")
    source_dynamic = root / candidate_freeze["source_artifact"]
    require(sha256(source_dynamic) == candidate_freeze["source_artifact_sha256"], "E05B candidate source changed")
    require(candidate_freeze["candidate_count"] == 15, "candidate count is not 15")
    require(candidate_freeze["pair_count"] == 10 and candidate_freeze["triple_count"] == 5, "candidate order counts changed")
    require(mode_freeze["locked_count"] == 15 and mode_freeze["no_lockable_count"] == 0, "not every candidate has a frozen mode")
    require(mode_freeze["holdout_inspected_before_freeze"] is False, "mode freeze admits holdout inspection")
    require(endpoint_freeze["coalition_outcomes_inspected_before_freeze"] is False, "endpoint freeze admits coalition inspection")
    limits_path = root / endpoint_freeze["limits_path"]
    trace_path = root / endpoint_freeze["trace_path"]
    require(sha256(limits_path) == endpoint_freeze["limits_sha256"], "stress limits changed after freeze")
    require(sha256(trace_path) == endpoint_freeze["trace_sha256"], "calibration trace changed after freeze")

    limits = pd.read_parquet(limits_path)
    tracks = pd.read_parquet(artifact / "holdout" / "e05c_mode_tracks.parquet")
    results = pd.read_parquet(artifact / "holdout" / "e05c_connected_pole_externalities.parquet")
    builds = pd.read_parquet(artifact / "holdout" / "e05c_failure_ledger.parquet")
    decisions = pd.read_parquet(artifact / "tables" / "e05c_candidate_decisions.parquet")
    require(len(limits) == 8 and limits["status"].eq("SUCCESS").all(), "baseline calibration is incomplete")
    require(limits["endpoint_reason"].eq("VOLTAGE_MIN").all(), "unexpected baseline endpoint reason")
    require(len(tracks) == 5760 and tracks["status"].eq("SUCCESS").all(), "mode-track table is incomplete")
    require(float(tracks["BMAC"].min()) >= 0.90, "a mode track violates the BMAC freeze")
    require(float(tracks["residual"].max()) <= 1e-10, "a mode track violates the residual freeze")
    require(len(results) == 1080 and results["status"].eq("SUCCESS").all(), "connected result population is incomplete")
    require(results["seed_id"].nunique() == 8 and results["tau"].nunique() == 9, "seed/tau grid mismatch")
    require(max(results["composition_closure_real"].abs().max(), results["composition_closure_imag"].abs().max()) <= 1e-10, "composition closure failed")
    require(len(builds) == 6984 and builds["status"].eq("SUCCESS").all(), "physical build ledger changed")
    require(int(results["stress_conditioned_material_delta_zeta"].sum()) == 0, "material damping points changed")
    require(int(results["hard_composition_failure"].sum()) == 0, "hard composition failures changed")
    require(int(results["screen_composition_failure"].sum()) == 0, "screen composition failures changed")
    require(abs(float(results["rho_comp"].max()) - 0.0021175800029241453) < 1e-14, "maximum rho_comp changed")
    require(len(decisions) == 15, "candidate decision table is incomplete")
    require(decisions["C3b_A"].eq("REJECTED").all(), "a C3b-A candidate unexpectedly passes")
    require(decisions["C3b_B"].eq("REJECTED").all(), "a C3b-B candidate unexpectedly passes")
    require(final["C3b_A_stress_conditioned_amplification"] == "REJECTED", "C3b-A final changed")
    require(final["C3b_B_material_margin_consequence"] == "REJECTED", "C3b-B final changed")
    require(final["C3a_nominal_material_damping_externality"] == "REJECTED_UNCHANGED", "C3a was reopened")
    require(final["e06_consideration_authorized"] is False, "E06 consideration was authorized")
    require(final["E06_executed"] is False, "E06 was executed")
    require(final["new_paraemt_scientific_runs"] == 0, "new ParaEMT science was recorded")
    e05 = json.loads((root / "artifacts/tx3/E05_mechanism_consequence/E05_FINAL_DECISION.json").read_text(encoding="utf-8"))
    e05b = json.loads((root / "artifacts/tx3/E05B_dynamic_spectral_shift/E05B_FINAL_DECISION.json").read_text(encoding="utf-8"))
    require(e05["C3a_mechanism_consequence"] == "REJECTED", "parent E05 decision changed")
    require(all(e05b[key] == "SUPPORTED" for key in ("D1_exact_schur_factorization", "D2_reproducible_finite_pole_externality", "D3_connected_contour_localization")), "parent E05B decisions changed")
    require(not (root / "artifacts/tx3/E06").exists(), "E06 artifact directory exists")
    for index, stem in enumerate((
        "stress_path_definition", "margin_conditioning", "all_candidate_heatmap", "composition_vs_actual",
        "damping_screen", "conditioning", "pair_vs_triple",
    ), start=1):
        for suffix in (".pdf", ".png"):
            path = artifact / "figures" / f"figure_E05C_{index}_{stem}{suffix}"
            require(path.is_file() and path.stat().st_size > 1000, f"missing/empty figure {path.name}")
    for name in REQUIRED_REPORTS:
        require((artifact / "reports" / name).is_file(), f"missing report {name}")
    versioned_claims = root / "reports/papers/tx3_connected_intervention_calculus/CLAIMS_EVIDENCE_MATRIX_E05C.md"
    require(versioned_claims.is_file(), "versioned E05C claim matrix missing")
    if versioned_claims.is_file():
        text = versioned_claims.read_text(encoding="utf-8")
        require("| C3b-A |" in text and "| C3b-B |" in text, "versioned claim matrix lacks C3b rows")
    package_manifest = artifact / "manifests" / "INDEPENDENT_REVIEW_PACKAGE_MANIFEST.json"
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
