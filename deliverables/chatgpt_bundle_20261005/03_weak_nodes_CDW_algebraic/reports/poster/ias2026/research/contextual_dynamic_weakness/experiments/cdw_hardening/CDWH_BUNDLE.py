# ruff: noqa: E501
"""H32: build CDW_HARDENING_REVIEW_BUNDLE_20260913.zip and CDW_PAPER_SUBMISSION_BUNDLE_20260913.zip.

Review bundle: hardening docs, theory, figures, results/hardening (summary tables, gates, ALT cases), the
top-level hardened CSVs, hardening and CDW code, logs, the paper folder, README_FOR_REVIEW_HARDENING.md and a
raw-store manifest (per-store file count, bytes and a digest of the sorted file hashes). Raw per-task
checkpoints, venvs, caches and git objects are excluded.

Submission bundle: the files needed to compile the manuscript and supplement, plus the PDFs.
"""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import csv
import hashlib
import json
import zipfile
from pathlib import Path

import _infra as I

DATE = "20260913"
PAPER = I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"
REVIEW = I.PROJECT / f"CDW_HARDENING_REVIEW_BUNDLE_{DATE}.zip"
SUBMIT = I.PROJECT / f"CDW_PAPER_SUBMISSION_BUNDLE_{DATE}.zip"
RAW_MAN = HI.RESULTS / "CDW_HARDENING_RAW_MANIFEST.csv"
SKIP_DIRS = {"__pycache__", ".git", ".venv"}
PAPER_JUNK = {".aux", ".log", ".out", ".blg", ".synctex.gz", ".toc"}


def raw_manifest():
    rows = []
    for d in sorted(I.RAW.glob("H_*")):
        files = sorted(d.glob("*.json"))
        h = hashlib.sha256()
        size = 0
        for p in files:
            b = p.read_bytes()
            size += len(b)
            h.update(p.name.encode())
            h.update(hashlib.sha256(b).digest())
        rows.append({"store": f"raw/{d.name}", "n_files": len(files), "bytes": size, "digest_sha256": h.hexdigest()})
    with RAW_MAN.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["store", "n_files", "bytes", "digest_sha256"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def _files(root: Path, pattern="*"):
    for p in sorted(root.rglob(pattern)):
        if p.is_file() and not (SKIP_DIRS & set(p.parts)):
            yield p


def review_bundle():
    n = 0
    with zipfile.ZipFile(REVIEW, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        def add(p: Path, arc: str):
            nonlocal n
            zf.write(p, arc)
            n += 1
        docs = [p for p in I.DOCS.glob("CDW_*") if p.is_file() and ("HARDENING" in p.name or p.name in (
            "CDW_NOVELTY_BOUNDARY.md", "CDW_REVIEWER_ATTACK.md", "CDW_REVIEWER1_FINAL.md", "CDW_REVIEWER2_FINAL.md", "CDW_AUTHOR_RESPONSE_TO_INTERNAL_REVIEW.md"))]
        for p in docs:
            add(p, f"docs/{p.name}")
        for p in _files(I.PROJECT / "theory"):
            add(p, str(p.relative_to(I.PROJECT)))
        for p in _files(HI.FIGS):
            add(p, str(p.relative_to(I.PROJECT)))
        for p in _files(HI.RESULTS):
            add(p, str(p.relative_to(I.PROJECT)))
        for name in ("CDW_HARDENED_CLAIM_MATRIX.csv", "CDW_HARDENING_MASTER_RESULTS.csv", "CDW_HARDENING_RUN_MANIFEST.csv", "CDW_LITERATURE_GAP_MATRIX.csv",
                     "CDW_E6_corridor_definitions.json"):
            p = I.RESULTS / name
            if p.exists():
                add(p, f"results/{name}")
        for sub in ("cdw_hardening", "cdw"):
            for p in _files(I.PROJECT / "experiments" / sub, "*.py"):
                add(p, str(p.relative_to(I.PROJECT)))
        for p in _files(HI.LOGS):
            add(p, str(p.relative_to(I.PROJECT)))
        for p in _files(PAPER):
            if p.suffix not in PAPER_JUNK:
                add(p, "paper/" + str(p.relative_to(PAPER)).replace("\\", "/"))
        readme = I.PROJECT / "README_FOR_REVIEW_HARDENING.md"
        if readme.exists():
            add(readme, "README_FOR_REVIEW_HARDENING.md")
    return n


def submission_bundle():
    n = 0
    keep = ["main.tex", "main.bbl", "main.pdf", "cdw_numbers.tex", "references.bib", "supplement.tex", "supplement.pdf", "README.md", "sec_corridors.tex"]
    keep += [p.name for p in PAPER.glob("tab_*.tex")]
    with zipfile.ZipFile(SUBMIT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for name in keep:
            p = PAPER / name
            if p.exists():
                zf.write(p, name)
                n += 1
        for p in sorted((PAPER / "figures").glob("*.pdf")):
            zf.write(p, f"figures/{p.name}")
            n += 1
    return n


def main():
    rows = raw_manifest()
    nr = review_bundle()
    ns = submission_bundle()
    out = {"review_bundle": str(REVIEW), "review_files": nr, "review_bytes": REVIEW.stat().st_size,
           "review_sha256": hashlib.sha256(REVIEW.read_bytes()).hexdigest(),
           "submission_bundle": str(SUBMIT), "submission_files": ns, "submission_bytes": SUBMIT.stat().st_size,
           "submission_sha256": hashlib.sha256(SUBMIT.read_bytes()).hexdigest(),
           "raw_stores": len(rows), "raw_files": sum(r["n_files"] for r in rows), "raw_bytes": sum(r["bytes"] for r in rows)}
    HI.write_json("CDWH_BUNDLES.json", out)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
