"""Assemble the final IAS2026 closure bundle and checksum manifest."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
ROOT = CAMPAIGN.parents[1]
bundle_name = f"IAS2026_BULLETPROOF_CLOSURE_{datetime.now().strftime('%Y%m%d_%H%M')}"
bundle_root = CAMPAIGN / "bundle" / bundle_name
bundle_root.mkdir(parents=True, exist_ok=False)

def copy_tree(name: str) -> None:
    source = CAMPAIGN / name
    if source.exists():
        shutil.copytree(source, bundle_root / name)

for directory in ["docs", "env", "src", "raw", "derived", "figures", "reports", "logs", "prereg", "git"]:
    copy_tree(directory)

(bundle_root / "paper").mkdir()
for filename in ["FINAL_EXECUTIVE_REPORT.pdf", "FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf"]:
    shutil.copy2(CAMPAIGN / "paper" / filename, bundle_root / "paper" / filename)

for filename in ["README.md"]:
    shutil.copy2(CAMPAIGN / filename, bundle_root / filename)

head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
(bundle_root / "git" / "final_head.txt").write_text(head + "\n", encoding="utf-8")

manifest = {
    "bundle": bundle_name,
    "created_local": datetime.now().isoformat(timespec="seconds"),
    "final_head": head,
    "base_tag": "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE",
    "branch": "research/ias2026-bulletproof-closure-v1",
    "entry_points": {
        "executive_report": "paper/FINAL_EXECUTIVE_REPORT.pdf",
        "paper": "paper/FINAL_PAPER.pdf",
        "poster": "paper/FINAL_POSTER_DRAFT.pdf",
        "audit": "reports/FINAL_BULLETPROOF_AUDIT.md",
        "dashboard": "reports/GATE_DASHBOARD.md",
        "claim_ledger": "docs/CLAIM_LEDGER.md",
        "reproduction": "src/python/run_all_python.py",
    },
    "labels": ["IEEE39_VALIDATED", "THEOREM_VALIDATED", "NONLINEAR_TDS_VALIDATED", "POWERDYNAMICS_VALIDATED", "RETROSPECTIVE", "REFUTED", "NOT_TESTED", "STOPPED_BY_GATE"],
    "scientific_disposition": "frozen nominal closure reproduced; stronger V9 transfer expectation refuted; dependent robustness and second-model claims not tested",
}
(bundle_root / "REPRODUCIBILITY_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

hash_lines = []
for path in sorted(p for p in bundle_root.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    hash_lines.append(f"{digest}  {path.relative_to(bundle_root).as_posix()}")
(bundle_root / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")

zip_path = CAMPAIGN / "bundle" / f"{bundle_name}.zip"
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(p for p in bundle_root.rglob("*") if p.is_file()):
        archive.write(path, path.relative_to(bundle_root.parent).as_posix())

print(zip_path)
print(f"files={len(hash_lines)}")
