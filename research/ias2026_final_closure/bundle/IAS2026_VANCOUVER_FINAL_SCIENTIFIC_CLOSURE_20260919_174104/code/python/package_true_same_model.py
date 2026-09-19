"""Package the corrected true-same-model/alternative-model/holdout audit."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

from campaign_root import campaign_root_from_argv


CAMPAIGN = campaign_root_from_argv()
REPO = CAMPAIGN.parents[1]
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
NAME = f"IAS2026_TRUE_SAME_MODEL_AND_HOLDOUT_{STAMP}"
BUNDLE = CAMPAIGN / "bundle" / NAME
BUNDLE.mkdir(parents=True, exist_ok=False)


def copy_file(relative: str, target: str | None = None) -> None:
    source = CAMPAIGN / relative
    if not source.is_file():
        return
    destination = BUNDLE / (target or relative)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(relative: str, target: str | None = None) -> None:
    source = CAMPAIGN / relative
    if not source.is_dir():
        return
    destination = BUNDLE / (target or relative)
    for path in source.rglob("*"):
        if not path.is_file() or path.suffix.lower() == ".pdf" or "__pycache__" in path.parts:
            continue
        destination_path = destination / path.relative_to(source)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination_path)


copy_tree("raw")
copy_tree("src", "code/src")
copy_tree("env/julia", "code/env/julia")
copy_tree("prereg")
for relative in (
    "README.md",
    "reports/GATE_DASHBOARD.md",
    "reports/GATE_DASHBOARD.csv",
    "reports/P1_GFL_PARITY.md",
    "reports/P2_POWERDYNAMICS_V4_STATUS.md",
    "reports/P3_COLLECTIVE_MECHANISM_STATUS.md",
    "reports/P4_JULIA_TDS_STATUS.md",
    "reports/P5_SECOND_GFL_STATUS.md",
    "reports/P6_BLIND_HOLDOUT_STATUS.md",
    "reports/P6_GENUINE_HOLDOUT_STATUS.md",
    "reports/P1_P6_EXECUTION_SUMMARY.md",
    "docs/CLAIM_LEDGER.md",
    "docs/TRUE_SAME_MODEL_RECONCILIATION.md",
    "docs/POWERDYNAMICS_RECONCILIATION.md",
    "derived/tables/TABLE_CLAIM_LEDGER.csv",
):
    copy_file(relative)

(BUNDLE / "REPRODUCE.md").write_text(
    "# Fresh-extraction reproduction\n\n"
    "This bundle intentionally excludes generated FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf. Run from the extracted bundle root.\n\n"
    "```powershell\n"
    "python code/src/python/run_p1_canonical.py --campaign-root .\n"
    "python code/src/python/run_p1_canonical_transfer.py --campaign-root .\n"
    "python code/src/python/compare_p1_transfer.py --campaign-root .\n"
    "python code/src/python/run_true_same_model_gate.py --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p1_gfl11.jl --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p2_pd39_portfolios.jl --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p5_simplegfldc_ieee39.jl --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p4_pd_tds.jl --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p6_genuine_holdout.jl --phase discovery --campaign-root .\n"
    "python code/src/python/prepare_p6_genuine_predictions.py --campaign-root .\n"
    "julia --project=code/env/julia code/src/julia/run_p6_genuine_holdout.jl --phase holdout --campaign-root .\n"
    "python code/src/python/score_p6_genuine_holdout.py --campaign-root .\n"
    "```\n\n"
    "P2 and P5 are alternative dynamic-model negative holdouts. The true same-model gate is explicitly STOPPED_BY_GATE until the exact custom SynchronousMachine/AVR/PSS model is ported and reconciled. P6 is reported as a single-class procedural holdout; no kappa or discrimination claim is made.\n",
    encoding="utf-8",
)

try:
    git_commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
except Exception:
    git_commit = "UNAVAILABLE"
manifest = {
    "bundle": NAME,
    "release_status": "NOT_FINAL_PAPER_OR_POSTER",
    "created_at_local": datetime.now().isoformat(timespec="seconds"),
    "git_commit_at_packaging": git_commit,
    "result_labels": {
        "P1": "CANONICAL_PYTHON_VS_JULIA_DEVICE_AND_DOCUMENTED_HARNESS_PASS",
        "P2": "FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT",
        "P3": "NO_BLOCKER_FOUND_IN_ALTERNATIVE_MODEL",
        "P4": "FRESH_OFFICIAL_POWERDYNAMICS_TDS",
        "P5": "FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT",
        "P6": "FRESH_ALTERNATIVE_DYNAMIC_MODEL_BLIND_HOLDOUT_SINGLE_CLASS",
        "G3": "STOPPED_BY_GATE",
    },
    "same_model_gate": "STOPPED_BY_GATE",
    "excluded": ["generated FINAL_PAPER.pdf", "generated FINAL_POSTER_DRAFT.pdf", "paper/poster rewrite"],
    "fresh_extraction_commands": "REPRODUCE.md",
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
