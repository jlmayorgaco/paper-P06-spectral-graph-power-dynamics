from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    artifact = root / "artifacts/tx3/critical_mode_baseline_audit"
    tables = artifact / "tables"
    protocol_path = artifact / "protocol/CRITICAL_MODE_BASELINE_AUDIT_FREEZE.json"

    def require(condition: bool, message: str) -> None:
        if not condition: errors.append(message)

    require(protocol_path.is_file(), "missing audit protocol freeze")
    if not protocol_path.is_file(): return errors
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    for relative, expected in protocol["source_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen audit source {relative}")
        if path.is_file(): require(sha256(path) == expected, f"frozen audit source changed: {relative}")
    for relative, expected in protocol["input_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen audit input {relative}")
        if path.is_file(): require(sha256(path) == expected, f"frozen audit input changed: {relative}")
    required = {
        "critical_modes.parquet", "critical_mode_state_anatomy.parquet", "critical_mode_subset_poles.parquet",
        "critical_mode_externalities.parquet", "audit_failure_ledger.parquet",
        "externality_vs_critical_mode_comparison.parquet", "critical_mode_coalition_summary.parquet",
    }
    for name in required: require((tables / name).is_file(), f"missing table {name}")
    if any(not (tables / name).is_file() for name in required): return errors
    modes = pd.read_parquet(tables / "critical_modes.parquet")
    anatomy = pd.read_parquet(tables / "critical_mode_state_anatomy.parquet")
    subsets = pd.read_parquet(tables / "critical_mode_subset_poles.parquet")
    results = pd.read_parquet(tables / "critical_mode_externalities.parquet")
    builds = pd.read_parquet(tables / "audit_failure_ledger.parquet")
    comparison = pd.read_parquet(tables / "externality_vs_critical_mode_comparison.parquet")
    summary = pd.read_parquet(tables / "critical_mode_coalition_summary.parquet")
    decision = json.loads((artifact / "CRITICAL_MODE_BASELINE_AUDIT_DECISION.json").read_text(encoding="utf-8"))
    require(len(modes) == 8 and modes["status"].eq("SUCCESS").all(), "critical mode population incomplete")
    require(modes["classification"].eq("synchronous_electromechanical").all(), "critical physical classification changed")
    require((modes["damping_ratio"] < 0.05).all(), "a baseline critical mode is not below 5%")
    require(len(anatomy) == 8 * 140, "state anatomy population incomplete")
    require(len(subsets) == 8 * 25 and subsets["status"].eq("SUCCESS").all(), "subset-track population incomplete")
    require(float(subsets["BMAC"].min()) >= 0.90, "BMAC tracking threshold violated")
    require(float(subsets["eigenpair_residual"].max()) <= 1e-10, "eigenpair residual threshold violated")
    require(len(results) == 8 * 15 and results["status"].eq("SUCCESS").all(), "connected-case population incomplete")
    require(max(float(results["composition_closure_real"].abs().max()), float(results["composition_closure_imag"].abs().max())) <= 1e-10, "complex composition closure failed")
    require(len(comparison) == 15 and len(summary) == 15, "comparison/summary population incomplete")
    require(set(builds["status"]) == {"SUCCESS"}, "physical/eigen build failure exists")
    require(decision["scope"] == "DESCRIPTIVE_POST_MORTEM_ONLY" and decision["claim_gate_created"] is False, "audit scope changed")
    require(decision["statuses_unchanged"] == {"C3a": "REJECTED", "C3b-A": "REJECTED", "C3b-B": "REJECTED", "C3c": "UNRESOLVED"}, "historical claim status changed")
    require(decision["mechanism_surgery_authorized"] is False and decision["E06_authorized"] is False, "forbidden mechanism/E06 authorization")
    require((artifact / "critical_mode_eigenvectors.npz").is_file(), "missing critical-mode eigenvectors")
    require((artifact / "TX3_CRITICAL_MODE_BASELINE_AUDIT.md").is_file(), "missing audit report")
    require((root / "reports/papers/tx3_connected_intervention_calculus/TX3_CRITICAL_MODE_BASELINE_AUDIT.md").is_file(), "missing manuscript-side audit report")
    stems = (
        "figure_1_critical_mode_participation_anatomy", "figure_2_critical_mode_frequency_damping",
        "figure_3_critical_mode_rho_all_coalitions", "figure_4_externality_vs_margin_mode",
    )
    for stem in stems:
        for suffix in (".pdf", ".png"):
            path = artifact / "figures" / f"{stem}{suffix}"
            require(path.is_file() and path.stat().st_size > 1000, f"missing/empty figure {path.name}")
    e05d = json.loads((root / "artifacts/tx3/E05D_weak_grid/E05D_FINAL_DECISION.json").read_text(encoding="utf-8"))
    require(e05d["terminal_status"] == "BASELINE_INELIGIBLE" and e05d["C3c_WEAK_GRID_MATERIAL_EXTERNALITY"] == "UNRESOLVED", "E05D terminal interpretation changed")
    require(not (root / "artifacts/tx3/E06").exists(), "E06 artifact directory exists")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(); errors = validate(args.root.resolve())
    print(json.dumps({"root": str(args.root.resolve()), "status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
