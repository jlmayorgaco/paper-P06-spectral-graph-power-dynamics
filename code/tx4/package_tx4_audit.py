"""Package corrected TX4 audit upload and reproducibility bundles."""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DELIVERABLES = ROOT / "deliverables"
REPORT_PDF = ROOT / "output" / "pdf" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.pdf"
FIGURES = ROOT / "results" / "TX4_AUDIT_FIGURES"


def add_file(zf: zipfile.ZipFile, path: Path, archive_name: str) -> None:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.name in {"FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf"}:
        raise RuntimeError(f"forbidden legacy PDF in bundle: {path}")
    zf.write(path, archive_name)


def write_bundle(path: Path, entries: list[tuple[Path, str]], bundle_kind: str) -> str:
    manifest = []
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for source, archive_name in entries:
            add_file(zf, source, archive_name)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            manifest.append({"archive_name": archive_name, "sha256": digest, "bytes": source.stat().st_size})
        manifest_payload = {
            "bundle_kind": bundle_kind,
            "branch": "research/tx4-robustness-statistics-audit",
            "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
            "no_push": True,
            "entries": manifest,
        }
        zf.writestr("MANIFEST.json", json.dumps(manifest_payload, indent=2) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    base_docs = [
        ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.md",
        ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_HANDOFF.md",
        ROOT / "docs" / "TX4_ROBUSTNESS_AUDIT_CURRENT_STATE.md",
        ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_PREREG.md",
        ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_DEVIATIONS.md",
        ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_EXECUTION_PLAN.md",
    ]
    result_names = [
        "TX4_ROBUSTNESS_HEADLINE_CORRECTED.json", "TX4_ROBUSTNESS_INVARIANT_CHECKS.csv",
        "TX4_QMC_SUMMARY.csv", "TX4_MC_SUMMARY.csv",
        "TX4_QMC_EXACT_BLOCKER_SUMMARY.csv", "TX4_MC_EXACT_BLOCKER_SUMMARY.csv",
        "TX4_QMC_DELTA_H4_TRUE_MINIMALITY.csv", "TX4_MC_DELTA_H4_TRUE_MINIMALITY.csv",
        "TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv", "TX4_EXACT_BLOCKER_CALIBRATION_BLOCKERS.csv",
        "TX4_AUDIT_H4_FLAG_COMPARISON.csv", "TX4_AUDIT_H4_FLAG_CONFUSION.csv",
        "TX4_AUDIT_COMPLEMENT_BUG.csv", "TX4_AUDIT_ETA_H4.csv", "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv",
        "TX4_ETA_LOCAL_COLLECTIVE_EXACT_SUMMARY.json", "TX4_SURROGATE_GATE.json",
        "TX4_SURROGATE_MINIMALITY_VALIDATION.csv", "TX4_SURROGATE_BLOCKER_ERRORS.csv",
        "TX4_SURROGATE_ALPHA_VALIDATION.csv", "TX4_ROBUSTNESS_DATA_PROVENANCE.csv",
    ]
    upload_entries: list[tuple[Path, str]] = [(REPORT_PDF, "report/TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.pdf")]
    upload_entries += [(path, f"docs/{path.name}") for path in base_docs]
    upload_entries += [(ROOT / "results" / name, f"results/{name}") for name in result_names]
    upload_entries += [(path, f"figures/{path.name}") for path in sorted(FIGURES.glob("F*.png"))]

    repro_entries = list(upload_entries)
    repro_entries += [
        (ROOT / "results" / "TX4_QMC_ALL16_EXACT.parquet", "results/TX4_QMC_ALL16_EXACT.parquet"),
        (ROOT / "results" / "TX4_MC_ALL16_EXACT.parquet", "results/TX4_MC_ALL16_EXACT.parquet"),
        (ROOT / "results" / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet", "results/TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet"),
        (ROOT / "results" / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz", "results/TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz"),
    ]
    scripts = [
        "tx4_audit_blockers.py", "audit_existing_robustness.py", "run_tx4_audit_exact_all16.py",
        "analyze_tx4_audit_exact.py", "run_tx4_eta_local_collective.py", "build_tx4_audit_outputs.py",
        "check_tx4_audit_invariants.py", "build_tx4_audit_pdf.py", "package_tx4_audit.py",
    ]
    repro_entries += [(ROOT / "code" / "tx4" / name, f"code/tx4/{name}") for name in scripts]
    repro_entries.append((ROOT / "code" / "tx4" / "tests" / "test_tx4_audit_blockers.py", "code/tx4/tests/test_tx4_audit_blockers.py"))

    upload = DELIVERABLES / "TX4_ROBUSTNESS_STATISTICS_AUDIT_CHATGPT_UPLOAD.zip"
    repro = DELIVERABLES / "TX4_ROBUSTNESS_STATISTICS_AUDIT_REPRO.zip"
    upload_sha = write_bundle(upload, upload_entries, "ChatGPT upload")
    repro_sha = write_bundle(repro, repro_entries, "reproducibility")
    (DELIVERABLES / "TX4_ROBUSTNESS_STATISTICS_AUDIT_SHA256SUMS.txt").write_text(
        f"{upload_sha}  {upload.name}\n{repro_sha}  {repro.name}\n", encoding="utf-8"
    )
    print(json.dumps({"upload": str(upload), "upload_sha256": upload_sha, "repro": str(repro), "repro_sha256": repro_sha, "upload_entries": len(upload_entries), "repro_entries": len(repro_entries)}, indent=2))


if __name__ == "__main__":
    main()

