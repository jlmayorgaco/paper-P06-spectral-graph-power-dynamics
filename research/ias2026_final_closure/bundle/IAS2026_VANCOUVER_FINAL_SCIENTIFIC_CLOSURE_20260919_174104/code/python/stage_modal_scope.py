"""Stage the corrected same-policy modal-scope audit into raw/modal_scope.

The committed campaign is replayable from raw/modal_scope. External staging is
an explicit operation used only to import the verified canonical audit once.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
DEST = CAMPAIGN / "raw" / "modal_scope"
FILES = [
    "TX4_FINAL_MODAL_SCOPE_V9_SUMMARY.json",
    "TX4_FINAL_MODAL_SCOPE_HEADLINE.json",
    "TX4_V9_GLOBAL_VS_EM_TRUTH.csv",
    "TX4_V9_EM_BAND_CONFUSION.csv",
    "TX4_V9_FALSE_SAFE_TAXONOMY.csv",
    "TX4_V9_SLOW_FALSE_SAFE_AUDIT.csv",
    "TX4_V9_BLIND_VS_FULL.csv",
    "TX4_V4_BLIND_VS_FULL.csv",
    "TX4_TRUE_LOCAL_VS_COLLECTIVE.csv",
    "TX4_TRUE_COLLECTIVE_SUBSET_FACTORS.csv",
    "TX4_Q_VS_RETURN_SIGN_AUDIT.csv",
    "TX4_PROPER_SUBSET_CLOSURE.csv",
    "TX4_15_PROPER_SUBSETS.csv",
    "TX4_MINIMALITY_SEPARATION.csv",
    "TX4_MODAL_MECHANISM.csv",
    "TX4_FINAL_MODAL_SCOPE_CLAIM_MATRIX.csv",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-mode", choices=("bundled", "external"), default="bundled")
    parser.add_argument("--source-root", type=Path, default=None)
    args = parser.parse_args()
    if args.source_mode == "bundled":
        source = DEST
    else:
        if args.source_root is None:
            raise SystemExit("--source-root is required with --source-mode external")
        source = args.source_root / "results"
    DEST.mkdir(parents=True, exist_ok=True)
    missing = [name for name in FILES if not (source / name).exists()]
    if missing:
        raise SystemExit("missing modal-scope sources:\n" + "\n".join(missing))
    if args.source_mode == "external":
        for name in FILES:
            shutil.copy2(source / name, DEST / name)
        report = args.source_root / "docs" / "TX4_FINAL_MODAL_SCOPE_REPORT.md"
        if report.exists():
            shutil.copy2(report, DEST / "TX4_FINAL_MODAL_SCOPE_REPORT.md")
    hashes = {}
    for path in sorted(DEST.iterdir()):
        if path.is_file():
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (DEST / "SOURCE_MANIFEST.json").write_text(json.dumps({
        "mode": args.source_mode,
        "source_root": str(args.source_root) if args.source_root else "bundled raw/modal_scope",
        "files": hashes,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"modal_scope_files={len(hashes)} mode={args.source_mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
