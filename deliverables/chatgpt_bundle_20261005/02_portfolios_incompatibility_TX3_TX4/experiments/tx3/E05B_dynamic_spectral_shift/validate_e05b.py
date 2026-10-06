from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    artifact = root / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    prereg = artifact / "preregistration"
    for name in (
        "E05B_GATE_FREEZE.json", "E05B_HOLDOUT_FREEZE.json", "E05B_HOLDOUT_MANIFEST.parquet",
        "E05B_NUMERICAL_FREEZE.json", "E05B_PREREGISTRATION.md", "E05B_CANDIDATE_CONTOUR_FREEZE.json",
    ):
        require((prereg / name).is_file(), f"missing preregistration/{name}")
    if errors:
        return errors
    numerical = json.loads((prereg / "E05B_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    freeze = json.loads((prereg / "E05B_CANDIDATE_CONTOUR_FREEZE.json").read_text(encoding="utf-8"))
    final = json.loads((artifact / "E05B_FINAL_DECISION.json").read_text(encoding="utf-8"))
    for relative, expected in numerical["source_hashes"].items():
        path = root / relative
        require(path.is_file(), f"missing frozen source {relative}")
        if path.is_file():
            require(sha256(path) == expected, f"frozen source hash mismatch: {relative}")
    development_path = artifact / "development" / "component_externalities.parquet"
    spectra_path = artifact / "spectra" / "development_spectra.npz"
    require(sha256(development_path) == freeze["development_component_sha256"], "development table changed after candidate freeze")
    require(sha256(spectra_path) == freeze["development_spectra_sha256"], "development spectra changed after candidate freeze")
    require(freeze["holdout_inspected_before_freeze"] is False, "candidate freeze admits prior holdout inspection")
    require(freeze["top_dynamic_pairs"] == ["A7-A8", "A3-A6"], "unexpected top dynamic pairs")
    require(freeze["top_dynamic_triples"] == ["A2-A7-A8", "A2-A5-A8"], "unexpected top dynamic triples")

    for population, op_count, vertex_count, case_count in (
        ("development", 24, 2232, 2016), ("holdout", 16, 1488, 1344)
    ):
        vertices = pd.read_parquet(artifact / population / "factor_vertices.parquet")
        cases = pd.read_parquet(artifact / population / "component_externalities.parquet")
        require(vertices["operating_point_id"].nunique() == op_count, f"{population} op count mismatch")
        require(len(vertices) == vertex_count, f"{population} vertex count mismatch")
        require(len(cases) == case_count, f"{population} case count mismatch")
        require(cases["evaluation_status"].eq("SUCCESS").all(), f"{population} has failed cases")
        require(cases["component_closure_residual"].abs().max() <= 1e-10, f"{population} component closure failed")
    development_vertices = pd.read_parquet(artifact / "development" / "factor_vertices.parquet")
    require(development_vertices["factorization_residual_fine"].abs().max() <= 1e-7, "integrated E03 factorization failed")
    audits = pd.concat(
        [pd.read_parquet(artifact / pop / "pointwise_factorization_audit.parquet") for pop in ("development", "holdout")],
        ignore_index=True,
    )
    require(len(audits) == 198, "pointwise audit row count mismatch")
    require(audits["residual"].abs().max() <= 1e-7, "pointwise factorization audit failed")
    holdout_freeze = json.loads((prereg / "E05B_HOLDOUT_FREEZE.json").read_text(encoding="utf-8"))
    require(len(holdout_freeze["operating_points"]) == 16, "holdout freeze does not contain 16 points")
    require(all(record["operating_point_id"].startswith("B") for record in holdout_freeze["operating_points"]), "holdout reuses E03 IDs")

    dynamic = pd.read_parquet(artifact / "tables" / "holdout_dynamic_reproducibility.parquet")
    require(int(((dynamic["order"] == 2) & dynamic["reproducible_dynamic_externality"]).sum()) == 10, "reproducible pair count mismatch")
    require(int(((dynamic["order"] == 3) & dynamic["reproducible_dynamic_externality"]).sum()) == 5, "reproducible triple count mismatch")
    contours = pd.read_parquet(artifact / "tables" / "contour_candidate_decisions.parquet")
    require(dict(zip(contours["coalition_key"], contours["status"], strict=True)) == {
        "A7-A8": "SUPPORTED", "A3-A6": "REJECTED", "A2-A7-A8": "SUPPORTED", "A2-A5-A8": "SUPPORTED"
    }, "contour decisions changed")
    require(final["D1_exact_schur_factorization"] == "SUPPORTED", "D1 not supported")
    require(final["D2_reproducible_finite_pole_externality"] == "SUPPORTED", "D2 not supported")
    require(final["D3_connected_contour_localization"] == "SUPPORTED", "D3 not supported")
    require(final["C3a"] == "REJECTED_UNCHANGED", "E05 C3a was reopened")
    require(final["material_modal_margin_externality"] == "NOT_ESTABLISHED", "materiality was overclaimed")
    require(final["progression_authorized"] is False, "progression was authorized")
    require(final["E05C_or_E06_executed"] is False, "E05C/E06 was executed")
    require(final["new_paraemt_scientific_runs"] == 0, "new ParaEMT science recorded")
    e05 = json.loads((root / "artifacts" / "tx3" / "E05_mechanism_consequence" / "E05_FINAL_DECISION.json").read_text(encoding="utf-8"))
    require(e05["C3a_mechanism_consequence"] == "REJECTED", "parent E05 negative decision changed")
    require(not (root / "artifacts" / "tx3" / "E06").exists(), "E06 artifact directory exists")
    for stem in (
        "figure_E05B_1_factor_anatomy", "figure_E05B_2_full_vs_components",
        "figure_E05B_3_dynamic_replication", "figure_E05B_4_contour_moments",
    ):
        for suffix in (".pdf", ".png"):
            path = artifact / "figures" / f"{stem}{suffix}"
            require(path.is_file() and path.stat().st_size > 1000, f"missing/empty figure {path.name}")
    for name in (
        "E05B_FACTOR_DECOMPOSITION.md", "E05B_DYNAMIC_HOLDOUT.md", "E05B_CONNECTED_SPECTRAL_SHIFT.md",
        "E05B_FAILURE_LEDGER.md", "E05B_DECISION.md", "E05B_CLAIMS_EVIDENCE_ADDENDUM.md",
        "INDEPENDENT_REVIEW_README.md",
    ):
        require((artifact / "reports" / name).is_file(), f"missing report {name}")
    claims = (artifact / "reports" / "E05B_CLAIMS_EVIDENCE_ADDENDUM.md").read_text(encoding="utf-8")
    require(all(f"| {claim} |" in claims for claim in ("C3b", "C3c", "C3d")), "E05B claims addendum lacks decision rows")
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
