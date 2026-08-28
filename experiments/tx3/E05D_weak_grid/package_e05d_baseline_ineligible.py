from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile
from datetime import UTC, datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
OUTPUT = ROOT / f"TX3_E05D_WEAK_GRID_REVIEW_CHATGPT_UPLOAD_{datetime.now(UTC).date().isoformat()}.zip"
MANIFEST = ARTIFACT / "manifests" / "PACKAGE_MANIFEST.json"
VALIDATION = ARTIFACT / "manifests" / "PACKAGE_SOURCE_VALIDATION.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def eligible(path: Path) -> bool:
    parts = set(path.relative_to(ROOT).parts)
    return path.is_file() and "__pycache__" not in parts and ".pytest_cache" not in parts and path.suffix not in {".pyc", ".pyo"} and path not in {OUTPUT, MANIFEST}


def collect() -> list[Path]:
    files: set[Path] = set()
    for relative in ("artifacts/tx3/E05D_weak_grid", "experiments/tx3/E05D_weak_grid", "dicgrid"):
        base = ROOT / relative
        files.update(path for path in base.rglob("*") if eligible(path))
    for relative in (
        "README.md", "LICENSE", "pyproject.toml", "configs/tx3_gfl_common.yaml",
        "systems/ieee39_tx3_gfl3/andes_case.xlsx",
        "experiments/tx3/E03_finite_externality/campaign_core.py",
        "experiments/tx3/E05C_margin_conditioned/margin_conditioned.py",
        "artifacts/tx3/E03_E04/preregistration/ACTION_SET_FREEZE.json",
        "artifacts/tx3/E05B_dynamic_spectral_shift/tables/holdout_dynamic_reproducibility.parquet",
        "artifacts/tx3/E05B_dynamic_spectral_shift/E05B_FINAL_DECISION.json",
        "artifacts/tx3/E05_mechanism_consequence/E05_FINAL_DECISION.json",
        "artifacts/tx3/E05C_margin_conditioned/E05C_FINAL_DECISION.json",
        "reports/papers/tx3_connected_intervention_calculus/CLAIMS_EVIDENCE_MATRIX.md",
        "reports/papers/tx3_connected_intervention_calculus/CLAIMS_EVIDENCE_MATRIX_E05C.md",
        "reports/papers/tx3_connected_intervention_calculus/CLAIMS_EVIDENCE_MATRIX_E05D.md",
        "tests/unit/test_e05d_weak_grid.py",
    ):
        path = ROOT / relative
        if not path.is_file(): raise FileNotFoundError(relative)
        files.add(path)
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def validate(root: Path) -> dict[str, object]:
    command = [sys.executable, str(root / "experiments/tx3/E05D_weak_grid/validate_e05d_baseline_ineligible.py"), "--root", str(root)]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    return {"return_code": completed.returncode, "stdout": completed.stdout, "stderr": completed.stderr, "status": "PASS" if completed.returncode == 0 else "FAIL"}


def main() -> int:
    source_validation = validate(ROOT)
    VALIDATION.write_text(json.dumps(source_validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if source_validation["status"] != "PASS": print(json.dumps(source_validation, indent=2)); return 1
    files = collect()
    records = [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in files]
    payload = {
        "package": OUTPUT.name, "created_utc": datetime.now(UTC).isoformat(),
        "scope": "TX3 E05D terminal baseline-ineligibility review; no coalition outcomes; stop before E06",
        "terminal_status": "BASELINE_INELIGIBLE", "file_count_excluding_manifest": len(records),
        "uncompressed_bytes_excluding_manifest": sum(record["bytes"] for record in records),
        "coalition_results_included": False, "new_paraemt_scientific_results_included": False,
        "E05E_results_included": False, "E06_results_included": False,
        "historical_package_references": {
            "E05B": "TX3_E05B_DYNAMIC_SPECTRAL_SHIFT_REVIEW_CHATGPT_UPLOAD_2026-08-27.zip",
            "E05C": "TX3_E05C_MARGIN_CONDITIONED_REVIEW_CHATGPT_UPLOAD_2026-08-27.zip",
        },
        "manifest_self_excluded_from_hash_list": True, "files": records,
    }
    MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with zipfile.ZipFile(OUTPUT, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files: archive.write(path, path.relative_to(ROOT).as_posix())
        archive.write(MANIFEST, MANIFEST.relative_to(ROOT).as_posix())
    with tempfile.TemporaryDirectory(prefix="tx3_e05d_terminal_review_") as temporary:
        extracted = Path(temporary)
        with zipfile.ZipFile(OUTPUT, "r") as archive: archive.extractall(extracted)
        extraction_validation = validate(extracted)
    if extraction_validation["status"] != "PASS": print(json.dumps(extraction_validation, indent=2)); return 1
    digest = sha256(OUTPUT)
    sidecar = OUTPUT.with_suffix(OUTPUT.suffix + ".sha256")
    sidecar.write_text(f"{digest}  {OUTPUT.name}\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS", "zip": str(OUTPUT), "zip_bytes": OUTPUT.stat().st_size,
        "zip_sha256": digest, "sidecar": str(sidecar), "file_count_excluding_manifest": len(records),
        "clean_extraction_validation": extraction_validation,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
