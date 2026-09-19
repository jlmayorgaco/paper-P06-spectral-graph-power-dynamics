"""Package the current corrective audit without mislabeling it as final.

The historical generated paper/poster PDFs are intentionally excluded. A
final release is a separate gated action after every required experiment has
passed or has been explicitly dispositioned.
"""
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
bundle_name = f"IAS2026_BULLETPROOF_CORRECTIVE_AUDIT_{datetime.now().strftime('%Y%m%d_%H%M')}"
bundle_root = CAMPAIGN / "bundle" / bundle_name
bundle_root.mkdir(parents=True, exist_ok=False)


def copy_tree(name: str) -> None:
    source = CAMPAIGN / name
    if source.exists():
        shutil.copytree(source, bundle_root / name)


for directory in ["docs", "env", "src", "raw", "derived", "figures", "reports", "logs", "prereg", "git"]:
    copy_tree(directory)

(bundle_root / "paper_source").mkdir()
shutil.copy2(CAMPAIGN / "src" / "python" / "build_pdfs.py", bundle_root / "paper_source" / "build_pdfs.py")
shutil.copy2(CAMPAIGN / "docs" / "CLAIM_LEDGER.md", bundle_root / "paper_source" / "CLAIM_LEDGER.md")

(bundle_root / "README_FIRST.md").write_text(
    "# Read first — IAS2026 corrective audit\n\n"
    "This is an interim corrective-audit bundle, not a final paper/poster release.\n"
    "The historical FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf are intentionally\n"
    "excluded and must not be used for Vancouver.\n\n"
    "Gate 0 and PowerDynamics Gate A pass in their declared scopes. Same-model\n"
    "GFL parity, robustness, second-model, new blind-holdout, and new Julia-TDS\n"
    "requirements remain incomplete. See FINAL_EXECUTIVE_REPORT.md,\n"
    "CLAIM_LEDGER.csv, and OPEN_LIMITATIONS.md.\n",
    encoding="utf-8",
)
(bundle_root / "FINAL_EXECUTIVE_REPORT.md").write_text(
    "# Executive report — corrective audit (not final)\n\n"
    "Release status: NOT_FINAL / STOPPED_BY_GATE for dependent final-release gates.\n\n"
    + (CAMPAIGN / "reports" / "FINAL_BULLETPROOF_AUDIT.md").read_text(encoding="utf-8"),
    encoding="utf-8",
)
shutil.copy2(CAMPAIGN / "derived" / "tables" / "TABLE_CLAIM_LEDGER.csv", bundle_root / "CLAIM_LEDGER.csv")
shutil.copy2(CAMPAIGN / "docs" / "OPEN_LIMITATIONS.md", bundle_root / "OPEN_LIMITATIONS.md")


def git_meta(repo: Path, rev: str) -> str:
    result = subprocess.run(["git", "-C", str(repo), "show", "-s", "--format=%H %cI %s", rev], capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else f"{rev} NOT_FOUND"


commit_lines = [
    "# Provenance commits (final release not authorized)",
    "prereg_freeze: " + git_meta(ROOT, "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE"),
    "prediction_freeze: " + git_meta(ROOT, "7c2a8fd60dde1bc9c77e65c2a360c50346d0bbdc"),
    "reveal: " + git_meta(ROOT, "96052a3e69964ff76a0fb47aa4767b2d8944f7ca"),
    "canonical_modal_scope_audit: " + git_meta(Path(r"C:\w\tx4modalaudit"), "9968de79"),
    "campaign_implementation: " + git_meta(ROOT, "f6db1e76"),
    "campaign_bundle_assembler: " + git_meta(ROOT, "c1e062fc"),
    "corrective_audit_current_head: " + git_meta(ROOT, "HEAD"),
]
(bundle_root / "git" / "commits_created.txt").write_text("\n".join(commit_lines) + "\n", encoding="utf-8")

head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
manifest = {
    "bundle": bundle_name,
    "release_status": "NOT_FINAL",
    "created_local": datetime.now().isoformat(timespec="seconds"),
    "current_head": head,
    "base_tag": "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE",
    "branch": "research/ias2026-bulletproof-closure-v1",
    "excluded": ["FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf", "FINAL_EXECUTIVE_REPORT.pdf"],
    "entry_points": {
        "read_first": "README_FIRST.md",
        "executive_report": "FINAL_EXECUTIVE_REPORT.md",
        "claim_ledger": "CLAIM_LEDGER.csv",
        "limitations": "OPEN_LIMITATIONS.md",
        "dashboard": "reports/GATE_DASHBOARD.md",
        "gate0": "reports/GATE0_REPORT.md",
        "powerdynamics_gate_a": "raw/powerdynamics/pd39_equilibrium_gate.md",
        "reproduction": "src/python/run_all_python.py",
    },
    "labels": ["PROVED", "NUMERICALLY_VERIFIED", "IEEE39_VALIDATED", "POWERDYNAMICS_VALIDATED", "SECOND_MODEL_VALIDATED", "NONLINEAR_TDS_VALIDATED", "SYNTHETIC_PILOT", "CONSTRUCTED_COUNTEREXAMPLE", "RETROSPECTIVE", "BLIND_HOLDOUT", "NOT_TESTED", "REFUTED", "UNRESOLVED", "STOPPED_BY_GATE"],
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
