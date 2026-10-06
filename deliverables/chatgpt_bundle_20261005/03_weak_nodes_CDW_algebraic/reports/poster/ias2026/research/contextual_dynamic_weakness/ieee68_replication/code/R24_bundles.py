# ruff: noqa: E501
"""R24 bundles. No virtualenvs, caches or git objects.

    python R24_bundles.py manifest   -> CDW68_RAW_MANIFEST.csv (every raw/ file: phase, bytes, sha256, task summary)
    python R24_bundles.py zip <HEAD> -> CDW68_REPLICATION_REVIEW_BUNDLE_20260913.zip and CDW68_RAW_RESULTS_20260913.zip
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import csv
import hashlib
import json
import sys
import zipfile

PAPER = R.RESEARCH.parents[2] / "papers" / "cdw_contextual_dynamic_weakness"
REPO = R.RESEARCH.parents[3]
MANIFEST = R.P68 / "CDW68_RAW_MANIFEST.csv"
SKIP = ("__pycache__", ".venv", ".git", ".pyc", ".zip")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def manifest():
    rows = []
    for p in sorted(R.RAW.rglob("*.json")):
        try:
            d = json.loads(p.read_text())
            t = d.get("task", {})
            summ = {k: t.get(k) for k in ("phase", "pid", "variant", "S", "draw", "target") if k in t}
            ok = d.get("ok")
        except Exception:  # noqa: BLE001
            summ, ok = {}, None
        rows.append({"path": p.relative_to(R.P68).as_posix(), "phase": p.parent.name, "bytes": p.stat().st_size, "sha256": sha(p), "ok": ok, "task": json.dumps(summ, sort_keys=True)})
    with open(MANIFEST, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", MANIFEST, len(rows), "files", sum(r["bytes"] for r in rows), "bytes")


def _add(z, p, arc, listing):
    if any(s in p.as_posix() for s in SKIP):
        return
    z.write(p, arc)
    listing.append((arc, p.stat().st_size, sha(p)))


def headline_raw():
    """Raw per-task records of the reference policy (P68_10 / B68_10), all 64 portfolios, both models."""
    out = []
    for p in sorted((R.RAW / "R06A").glob("*.json")):
        t = json.loads(p.read_text())["task"]
        if t.get("pid") == "P68_10" and t.get("variant") == "REAL":
            out.append(p)
    for p in sorted((R.RAW / "B68C").glob("*.json")):
        t = json.loads(p.read_text()).get("task", {})
        if t.get("pid", t.get("cond")) in ("B68_10",) or str(t.get("cid", "")) == "B68_10":
            out.append(p)
    return out


def zips(head):
    listing = []
    rb = R.P68 / "CDW68_REPLICATION_REVIEW_BUNDLE_20260913.zip"
    with zipfile.ZipFile(rb, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for sub in ("docs", "inputs", "code", "results", "figures", "logs"):
            for p in sorted((R.P68 / sub).rglob("*")):
                if p.is_file():
                    _add(z, p, f"ieee68_replication/{p.relative_to(R.P68).as_posix()}", listing)
        for f in (".gitignore", "CDW68_RAW_MANIFEST.csv"):
            _add(z, R.P68 / f, f"ieee68_replication/{f}", listing)
        hr = headline_raw()
        for p in hr:
            _add(z, p, f"ieee68_replication/{p.relative_to(R.P68).as_posix()}", listing)
        for p in sorted(PAPER.rglob("*")):
            if p.is_file() and p.suffix in (".tex", ".bib", ".bbl", ".pdf", ".md"):
                _add(z, p, f"paper/{p.relative_to(PAPER).as_posix()}", listing)
        readme = (f"CDW68 replication review bundle (2026-09-13)\nbranch research/cdw-ieee68-replication, HEAD {head}\n"
                  "decision: CASE B (REV68A FAIL, REV68B FAIL, RANK68A PASS, RANK68B PASS)\n"
                  "start with ieee68_replication/docs/CDW68_FINAL_REPLICATION_REPORT.pdf; paper in paper/main.pdf\n"
                  f"raw headline records: {len(hr)} per-task JSON files of the reference policy (P68_10, B68_10); full raw data in CDW68_RAW_RESULTS_20260913.zip\n"
                  "no virtualenvs, caches or git objects\n")
        z.writestr("README.txt", readme)
        z.writestr("BUNDLE_MANIFEST.csv", "path,bytes,sha256\n" + "\n".join(f"{a},{b},{c}" for a, b, c in listing) + "\n")
    print("wrote", rb, len(listing), "files", rb.stat().st_size, "bytes; headline raw", len(hr))
    rz = R.P68 / "CDW68_RAW_RESULTS_20260913.zip"
    n = 0
    with zipfile.ZipFile(rz, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in sorted(R.RAW.rglob("*.json")):
            z.write(p, f"ieee68_replication/{p.relative_to(R.P68).as_posix()}")
            n += 1
        z.write(MANIFEST, "ieee68_replication/CDW68_RAW_MANIFEST.csv")
    print("wrote", rz, n, "files", rz.stat().st_size, "bytes")


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        manifest()
    else:
        zips(sys.argv[2])
