"""Create upload and reproducibility bundles for the TX4 campaign."""

from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DELIVERABLES = ROOT / "deliverables"
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
DELIVERABLES.mkdir(parents=True, exist_ok=True)


def add(zf: zipfile.ZipFile, path: Path, prefix: str) -> None:
    zf.write(path, f"{prefix}/{path.relative_to(ROOT).as_posix()}")


upload_files = [
    DOCS / "TX4_ROBUSTNESS_FINAL_REPORT.pdf",
    DOCS / "TX4_ROBUSTNESS_FINAL_REPORT.md",
    DOCS / "TX4_ROBUSTNESS_STATISTICAL_APPENDIX.pdf",
    DOCS / "TX4_ROBUSTNESS_STATISTICAL_APPENDIX.md",
    DOCS / "TX4_IAS_ROBUSTNESS_POSTER_SUMMARY.md",
    DOCS / "TX4_ROBUSTNESS_CHATGPT_HANDOFF.md",
    RESULTS / "TX4_ROBUSTNESS_HEADLINE.json",
    RESULTS / "TX4_ROBUSTNESS_CLAIM_MATRIX.csv",
    RESULTS / "TX4_QMC_SUMMARY.csv",
    RESULTS / "TX4_MC_SUMMARY.csv",
    RESULTS / "TX4_METHOD_COMPARISON.csv",
    RESULTS / "TX4_RETURN_ROBUSTNESS.csv",
    RESULTS / "TX4_MODE_ROBUSTNESS.csv",
]

repro_files = upload_files + [
    DOCS / "TX4_ROBUSTNESS_CURRENT_STATE.md",
    DOCS / "TX4_ROBUSTNESS_PREREG.md",
    DOCS / "TX4_ROBUSTNESS_STATISTICAL_PLAN.md",
    DOCS / "TX4_ROBUSTNESS_CLAIMS_V0.csv",
    DOCS / "TX4_ROBUSTNESS_DEVIATIONS.md",
    DOCS / "TX4_ROBUSTNESS_EXECUTION_PLAN.md",
    RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz",
    RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet",
    RESULTS / "TX4_ROBUSTNESS_MASTER_METADATA.json",
    RESULTS / "TX4_ROBUSTNESS_EXACT_CALIBRATION.csv",
    RESULTS / "TX4_ROBUSTNESS_H4_CHECKPOINTS.csv",
    RESULTS / "TX4_1D_SWEEPS.csv",
    RESULTS / "TX4_2D_PHASE_MAPS.parquet",
    RESULTS / "TX4_MORRIS.csv",
    RESULTS / "TX4_SOBOL.csv",
    RESULTS / "TX4_LOGISTIC_MODEL.csv",
    RESULTS / "TX4_GSTAR_DISTRIBUTION.csv",
    RESULTS / "TX4_TDS_ROBUSTNESS.csv",
    RESULTS / "TX4_JULIA_CROSSCODE_RANDOM.csv",
    ROOT / "code" / "tx4" / "run_tx4_robustness.py",
    ROOT / "code" / "tx4" / "run_tx4_robustness_campaign.py",
    ROOT / "code" / "tx4" / "run_tx4_h4_checkpoint.py",
    ROOT / "code" / "tx4" / "assemble_tx4_robustness.py",
    ROOT / "code" / "tx4" / "normalize_tx4_transverse_endpoints.py",
    ROOT / "code" / "tx4" / "analyze_tx4_robustness.py",
    ROOT / "code" / "tx4" / "build_tx4_robustness_reports.py",
    ROOT / "code" / "tx4" / "package_tx4_robustness.py",
]

for path in upload_files + repro_files:
    if not path.exists():
        raise FileNotFoundError(path)

upload = DELIVERABLES / "TX4_ROBUSTNESS_CHATGPT_UPLOAD.zip"
repro = DELIVERABLES / "TX4_ROBUSTNESS_REPRO.zip"
for target in (upload, repro):
    if target.exists():
        target.unlink()

with zipfile.ZipFile(upload, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
    for path in upload_files:
        add(zf, path, "TX4_ROBUSTNESS_CHATGPT_UPLOAD")

with zipfile.ZipFile(repro, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
    for path in repro_files:
        add(zf, path, "TX4_ROBUSTNESS_REPRO")

print(upload)
print(repro)
print("upload_sha256", hashlib.sha256(upload.read_bytes()).hexdigest())
print("repro_sha256", hashlib.sha256(repro.read_bytes()).hexdigest())
