"""Package the two requested TX4 deliverables without repository metadata."""

from __future__ import annotations

import zipfile
from pathlib import Path

from campaign_root import campaign_root_from_argv

ROOT = campaign_root_from_argv()
DELIVERABLES = ROOT / "deliverables"


def add_file(zf: zipfile.ZipFile, path: Path, prefix: str) -> None:
    if path.is_file():
        zf.write(path, f"{prefix}/{path.relative_to(ROOT).as_posix()}")


def add_glob(zf: zipfile.ZipFile, pattern: str, prefix: str) -> None:
    for path in sorted(ROOT.glob(pattern)):
        if not path.is_file():
            continue
        if any(part in {".git", "__pycache__", ".julia", "depot", "tmp"} for part in path.parts):
            continue
        if path.suffix in {".pyc", ".ji"}:
            continue
        add_file(zf, path, prefix)


upload = DELIVERABLES / "TX4_EXACT_P4_JULIA_CHATGPT_UPLOAD.zip"
with zipfile.ZipFile(upload, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    prefix = "TX4_EXACT_P4_JULIA_CHATGPT_UPLOAD"
    for path in [
        DELIVERABLES / "TX4_EXACT_P4_JULIA_FINAL_REPORT.pdf",
        ROOT / "reports" / "TX4_EXACT_P4_JULIA_HANDOFF.md",
        ROOT / "results" / "TX4_EXACT_P4_HEADLINE.json",
        ROOT / "results" / "TX4_EXACT_P4_CLAIMS.csv",
        ROOT / "results" / "TX4_EXACT_P4_V4_CROSSCODE.csv",
        ROOT / "results" / "TX4_EXACT_P4_MODE_MATCH.csv",
        ROOT / "results" / "TX4_EXACT_P4_STATE_COUNTS.csv",
        ROOT / "docs" / "TX4_GFL11_PYTHON_VS_JULIA_EQUATION_AUDIT.md",
        ROOT / "docs" / "TX4_P4_SG_POLICY_AUDIT.md",
        ROOT / "figures" / "TX4_F1_ALPHA_CROSSCODE.png",
        ROOT / "figures" / "TX4_F2_FREQUENCY_CROSSCODE.png",
        ROOT / "figures" / "TX4_F3_H4_VOLTAGE_MODE.png",
        ROOT / "figures" / "TX4_F4_STATE_COUNTS.png",
    ]:
        add_file(zf, path, prefix)

repro = DELIVERABLES / "TX4_EXACT_P4_JULIA_REPRO.zip"
with zipfile.ZipFile(repro, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    prefix = "TX4_EXACT_P4_JULIA_REPRO"
    for pattern in [
        "code/tx4/*.py",
        "code/tx4/*.jl",
        "docs/TX4_EXACT_P4_*.md",
        "docs/TX4_EXACT_P4_*.csv",
        "docs/TX4_GFL11_PYTHON_VS_JULIA_EQUATION_AUDIT.md",
        "docs/TX4_P4_SG_POLICY_AUDIT.md",
        "reports/TX4_EXACT_P4_*.md",
        "results/TX4_*.csv",
        "results/TX4_*.json",
        "raw/reconciliation/*.csv",
        "raw/tx4_exact_p4_julia/*.csv",
        "figures/TX4_*.png",
        "research/ias2026_last_validation/code/python/ibr_cycles/**/*.py",
        "research/ias2026_last_validation/raw/true_same_model/canonical_source/**/*.py",
        "research/ias2026_last_validation/raw/true_same_model/canonical_source/**/*.json",
        "deliverables/TX4_EXACT_P4_JULIA_FINAL_REPORT.tex",
    ]:
        add_glob(zf, pattern, prefix)
    add_file(zf, DELIVERABLES / "TX4_EXACT_P4_JULIA_FINAL_REPORT.pdf", prefix)

print(upload)
print(repro)
