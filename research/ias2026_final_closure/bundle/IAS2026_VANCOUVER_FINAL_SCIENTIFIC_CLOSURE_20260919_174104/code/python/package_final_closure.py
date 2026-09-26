"""Package the final closure evidence without legacy paper/poster PDFs."""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
NAME = f"IAS2026_VANCOUVER_FINAL_SCIENTIFIC_CLOSURE_{STAMP}"
BUNDLE = ROOT / "bundle" / NAME
BUNDLE.mkdir(parents=True, exist_ok=False)


def copy_file(relative: str, target: str | None = None) -> None:
    source = ROOT / relative
    if not source.is_file():
        return
    destination = BUNDLE / (target or relative)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(relative: str) -> None:
    source = ROOT / relative
    if not source.is_dir():
        return
    for path in source.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        if path.name in {"FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf"}:
            raise RuntimeError(f"forbidden legacy PDF in package source: {path}")
        destination = BUNDLE / path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)


for directory in ("raw", "code", "env", "prereg", "reports", "docs", "derived", "figures", "logs"):
    copy_tree(directory)
copy_file("README.md")
copy_file("reports/FINAL_EXECUTIVE_REPORT.pdf", "FINAL_EXECUTIVE_REPORT.pdf")

reproduction = """# Fresh extraction reproduction

This package excludes FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf. Run from the extracted bundle root.

```powershell
python code/python/run_p1_canonical.py --campaign-root .
python code/python/run_p1_canonical_transfer.py --campaign-root .
python code/python/compare_p1_transfer.py --campaign-root .
python code/python/run_true_same_model_gate.py --campaign-root .
julia --project=env/julia code/julia/run_p1_gfl11.jl --campaign-root .
python code/python/build_executive_report_pdf.py --campaign-root .
```

The true same-model gate is intentionally STOPPED_BY_GATE until the exact custom synchronous-machine/AVR/PSS Julia port is completed. Alternative-model results remain model-conditioned negative evidence.
"""
(BUNDLE / "REPRODUCE.md").write_text(reproduction, encoding="utf-8")

manifest = {
    "bundle": NAME,
    "created_at_local": datetime.now().isoformat(timespec="seconds"),
    "release_status": "CLOSURE_PACKAGE_WITH_STOPPED_GATES",
    "git_commit_at_packaging": "ef59780d39ef87cb53ecd0a872ac92e469d79748",
    "branch": "research/ias2026-final-award-closure",
    "statuses": {
        "G0": "FROZEN_ARTIFACT_REPRODUCED",
        "G1": "NUMERICALLY_VERIFIED",
        "G2_G4": "STOPPED_BY_GATE",
        "G5": "NOT_TESTED",
        "G6": "NUMERICALLY_VERIFIED",
        "G7": "SECOND_MODEL_VALIDATED",
        "G8": "BLIND_HOLDOUT",
        "G9_G12": "NOT_TESTED",
        "G13_G15": "RETROSPECTIVE",
        "G16": "NOT_TESTED",
    },
    "excluded": ["FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf", "paper/poster rewrite"],
    "fresh_extraction_commands": "REPRODUCE.md",
}
(BUNDLE / "Project_Manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

for path in BUNDLE.rglob("*"):
    if path.is_file() and path.name in {"FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf"}:
        raise RuntimeError(f"forbidden PDF entered package: {path}")

hash_lines = []
for path in sorted(path for path in BUNDLE.rglob("*") if path.is_file()):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    hash_lines.append(f"{digest}  {path.relative_to(BUNDLE).as_posix()}")
(BUNDLE / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")

zip_path = ROOT / "bundle" / f"{NAME}.zip"
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(path for path in BUNDLE.rglob("*") if path.is_file()):
        archive.write(path, path.relative_to(BUNDLE.parent).as_posix())
print(zip_path)
