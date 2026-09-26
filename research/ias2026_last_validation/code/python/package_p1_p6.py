"""Package the executed P1-P6 audit without paper/poster release artifacts."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[2]
REPO = CAMPAIGN.parents[1]
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
NAME = f"IAS2026_P1_P6_EXECUTION_{STAMP}"
BUNDLE = CAMPAIGN / "bundle" / NAME
BUNDLE.mkdir(parents=True, exist_ok=False)


def copy_file(relative: str, target: str | None = None) -> None:
    source = CAMPAIGN / relative
    destination = BUNDLE / (target or relative)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(relative: str, target: str | None = None) -> None:
    source = CAMPAIGN / relative
    destination = BUNDLE / (target or relative)
    shutil.copytree(source, destination, dirs_exist_ok=True)


copy_tree("raw/gfl11")
for section in ("p2", "p3", "p4", "p5", "p6"):
    copy_tree(f"raw/{section}")
copy_tree("src", "code/src")
copy_tree("env/julia", "code/env/julia")
copy_tree("prereg")
for relative in (
    "reports/GATE_DASHBOARD.md",
    "reports/P1_GFL_PARITY.md",
    "reports/P2_POWERDYNAMICS_V4_STATUS.md",
    "reports/P3_COLLECTIVE_MECHANISM_STATUS.md",
    "reports/P4_JULIA_TDS_STATUS.md",
    "reports/P5_SECOND_GFL_STATUS.md",
    "reports/P6_BLIND_HOLDOUT_STATUS.md",
    "reports/P1_P6_EXECUTION_SUMMARY.md",
    "reports/OPEN_LIMITATIONS.md",
    "docs/CLAIM_LEDGER.md",
    "docs/POWERDYNAMICS_RECONCILIATION.md",
    "derived/tables/TABLE_CLAIM_LEDGER.csv",
):
    if (CAMPAIGN / relative).is_file():
        copy_file(relative)

git_commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
manifest = {
    "bundle": NAME,
    "release_status": "NOT_FINAL_PAPER_OR_POSTER",
    "created_at_local": datetime.now().isoformat(timespec="seconds"),
    "git_commit_at_packaging": git_commit,
    "gates": {
        "P1": "PASS_CANONICAL_PYTHON_JULIA_DEVICE_AND_NETWORK",
        "P2": "PASS_16_OF_16",
        "P3": "NO_BLOCKER_FOUND",
        "P4": "NO_BLOCKER_FOUND",
        "P5": "PASS_SIMPLEGFLDC",
        "P6": "PASS_BLIND_HOLDOUT",
    },
    "required_evidence": [
        "canonical source/spec snapshot and SHA manifest",
        "corrected PowerDynamics current-source harness",
        "canonical and Julia transfer CSVs plus comparison",
        "P2 census/Hasse/mode shapes",
        "P3/P4 negative dependency results",
        "P5 source/config/result",
        "P6 prereg hash, committed predictions, reveal metrics",
        "dashboard, negative results, claim ledger, code, Project/Manifest",
    ],
    "excluded": ["generated FINAL_PAPER.pdf", "generated FINAL_POSTER_DRAFT.pdf", "paper/poster rewrite"],
}
(BUNDLE / "Project_Manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

hash_lines = []
for path in sorted(p for p in BUNDLE.rglob("*") if p.is_file()):
    hash_lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(BUNDLE).as_posix()}")
(BUNDLE / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
zip_path = CAMPAIGN / "bundle" / f"{NAME}.zip"
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(p for p in BUNDLE.rglob("*") if p.is_file()):
        archive.write(path, path.relative_to(BUNDLE.parent).as_posix())
print(zip_path)
