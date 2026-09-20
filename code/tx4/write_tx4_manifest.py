"""Write a compact reproducibility manifest for the exact TX4 bundle."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from campaign_root import campaign_root_from_argv

ROOT = campaign_root_from_argv()

patterns = [
    "code/tx4/*.py",
    "code/tx4/*.jl",
    "docs/TX4_EXACT_P4_*.md",
    "docs/TX4_EXACT_P4_*.csv",
    "docs/TX4_GFL11_PYTHON_VS_JULIA_EQUATION_AUDIT.md",
    "docs/TX4_P4_SG_POLICY_AUDIT.md",
    "reports/TX4_EXACT_P4_*.md",
    "results/TX4_*.csv",
    "results/TX4_*.json",
    "figures/TX4_*.png",
    "deliverables/TX4_EXACT_P4_JULIA_FINAL_REPORT.tex",
    "deliverables/TX4_EXACT_P4_JULIA_FINAL_REPORT.pdf",
]
paths: set[Path] = set()
for pattern in patterns:
    paths.update(ROOT.glob(pattern))

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

files = [
    {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": digest(path)}
    for path in sorted(paths)
    if path.is_file()
]
git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
git_branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
python_version = sys.version.split()[0]
try:
    julia_version = subprocess.check_output(["julia", "--version"], text=True).strip()
except Exception:
    julia_version = "unavailable"

manifest = {
    "bundle": "TX4_EXACT_P4_JULIA_REPRO",
    "git_branch": git_branch,
    "git_commit_at_manifest": git_commit,
    "python_version": python_version,
    "julia_version": julia_version,
    "policy": {"g": 0.03625, "k": 1.425, "t": 1.5, "h": 1.0},
    "gfl_states": 11,
    "target_buses": [30, 33, 35, 37],
    "files": files,
}
(ROOT / "results" / "TX4_EXACT_P4_REPRO_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"git_commit": git_commit, "files": len(files)}, indent=2))
