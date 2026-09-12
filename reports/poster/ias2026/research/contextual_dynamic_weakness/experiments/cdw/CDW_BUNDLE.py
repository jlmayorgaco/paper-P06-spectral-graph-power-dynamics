# ruff: noqa: E501
"""Phase 29: build the review bundle(s).

A. Compact review bundle: docs, theory, figures, results (summary tables), scripts,
   manifests, logs, README_FOR_REVIEW.md. Excludes venvs, git objects, caches, raw/.
B. Raw data bundle: raw/ (per-task checkpoints) compressed per phase, or a manifest if
   too large for a single archive.
"""

from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from pathlib import Path

import _infra as I

DATE = "20260912"
COMPACT = I.PROJECT / f"CDW_REVIEW_BUNDLE_{DATE}.zip"
RAW_MANIFEST = I.RESULTS / "CDW_RAW_RESULTS_MANIFEST.csv"
RAW_DIR = I.PROJECT / f"CDW_RAW_RESULTS_{DATE}"


def _add(zf: zipfile.ZipFile, path: Path, arc: str):
    zf.write(path, arc)


def compact_bundle():
    with zipfile.ZipFile(COMPACT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for sub in ("docs", "theory", "figures"):
            d = I.PROJECT / sub
            if not d.exists():
                continue
            for p in d.rglob("*"):
                if p.is_file():
                    _add(zf, p, str(p.relative_to(I.PROJECT)))
        for p in I.RESULTS.rglob("*"):
            if p.is_file() and p.suffix in (".csv", ".json"):
                _add(zf, p, str(p.relative_to(I.PROJECT)))
        for p in (I.PROJECT / "experiments").rglob("*.py"):
            if "__pycache__" not in p.parts:
                _add(zf, p, str(p.relative_to(I.PROJECT)))
        for p in I.LOGS.glob("*.log"):
            _add(zf, p, str(p.relative_to(I.PROJECT)))
        readme = I.PROJECT / "README_FOR_REVIEW.md"
        if readme.exists():
            _add(zf, readme, "README_FOR_REVIEW.md")
    return COMPACT


def raw_manifest_and_archives():
    rows = []
    RAW_DIR.mkdir(exist_ok=True)
    for phase_dir in sorted(I.RAW.iterdir()):
        if not phase_dir.is_dir():
            continue
        files = sorted(phase_dir.glob("*.json"))
        total_bytes = sum(f.stat().st_size for f in files)
        arc_path = RAW_DIR / f"CDW_RAW_{phase_dir.name}.zip"
        with zipfile.ZipFile(arc_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for f in files:
                zf.write(f, f"{phase_dir.name}/{f.name}")
        rows.append({"phase": phase_dir.name, "n_files": len(files), "raw_bytes": total_bytes,
                     "archive": str(arc_path.relative_to(I.PROJECT)), "archive_bytes": arc_path.stat().st_size,
                     "archive_sha256": hashlib.sha256(arc_path.read_bytes()).hexdigest()})
    with RAW_MANIFEST.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["phase", "n_files", "raw_bytes", "archive", "archive_bytes", "archive_sha256"])
        for r in rows:
            w.writerow([r["phase"], r["n_files"], r["raw_bytes"], r["archive"], r["archive_bytes"], r["archive_sha256"]])
    return rows


def main():
    rows = raw_manifest_and_archives()
    path = compact_bundle()
    print(json.dumps({"compact_bundle": str(path), "compact_bytes": path.stat().st_size,
                       "raw_phases": len(rows), "raw_total_bytes": sum(r["raw_bytes"] for r in rows)}, indent=1))


if __name__ == "__main__":
    main()
