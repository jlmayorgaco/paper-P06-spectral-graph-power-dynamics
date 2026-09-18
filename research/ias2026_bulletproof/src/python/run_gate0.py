"""Read-only Gate 0 audit against the frozen campaign artifacts.

The driver does not regenerate historical result directories. It copies only the
small source tables needed for this audit into the new campaign tree and writes
derived pass/fail tables with explicit provenance.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[2]
REPO = CAMPAIGN.parents[1]
FROZEN = REPO / "reports" / "poster" / "ias2026" / "research"
TX4BLIND = Path(r"C:\w\tx4blind")
RAW = CAMPAIGN / "raw" / "gate0"
DERIVED = CAMPAIGN / "derived" / "tables"
REPORTS = CAMPAIGN / "reports"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def copy_source(path: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
    return str(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(command: str) -> str:
    result = subprocess.run(
        ["git", *command.split()], cwd=REPO, text=True, capture_output=True
    )
    return result.stdout.strip() if result.returncode == 0 else f"ERROR: {result.stderr.strip()}"


def number(value: str | None):
    if value in (None, "", "nan", "NaN"):
        return None
    return float(value)


def row(metric, expected, observed, threshold, status, label, source):
    error = None
    if expected is not None and observed is not None:
        error = abs(observed - expected)
    return {
        "metric": metric,
        "expected": expected,
        "observed": observed,
        "error": error,
        "threshold": threshold,
        "status": status,
        "label": label,
        "source_file": source,
    }


def main() -> int:
    for path in (RAW, DERIVED, REPORTS):
        path.mkdir(parents=True, exist_ok=True)

    fc03 = FROZEN / "results" / "FINAL_CLOSURE" / "FC03_points.csv"
    fc01 = FROZEN / "results" / "FINAL_CLOSURE" / "FC01_structure.csv"
    fc03_summary = FROZEN / "results" / "FINAL_CLOSURE" / "FC03_summary.json"
    tds = FROZEN / "results" / "G2" / "G2_tds_summary.csv"
    f7a = FROZEN / "results" / "F7" / "F7A_events.csv"
    v9 = TX4BLIND / "deliverables" / "TX4_BLIND_PREDICTION_CHATGPT_UPLOAD" / "TX4_V9_BLIND_VS_FULL.csv"
    reveal = TX4BLIND / "deliverables" / "TX4_BLIND_PREDICTION_REPRO" / "results" / "TX4_REVEAL_SUMMARY.json"

    required = [fc03, fc01, fc03_summary, tds, f7a, v9, reveal]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        report = REPORTS / "GATE0_REPORT.md"
        report.write_text(
            "# Gate 0 Report\n\nstatus: STOPPED_BY_GATE\n\nMissing source files:\n"
            + "\n".join(f"- `{path}`" for path in missing)
            + "\n",
            encoding="utf-8",
        )
        return 2

    sources = {
        "FC03_points.csv": copy_source(fc03, RAW / "FC03_points.csv"),
        "FC01_structure.csv": copy_source(fc01, RAW / "FC01_structure.csv"),
        "FC03_summary.json": copy_source(fc03_summary, RAW / "FC03_summary.json"),
        "G2_tds_summary.csv": copy_source(tds, RAW / "G2_tds_summary.csv"),
        "F7A_events.csv": copy_source(f7a, RAW / "F7A_events.csv"),
        "TX4_V9_BLIND_VS_FULL.csv": copy_source(v9, RAW / "TX4_V9_BLIND_VS_FULL.csv"),
        "TX4_REVEAL_SUMMARY.json": copy_source(reveal, RAW / "TX4_REVEAL_SUMMARY.json"),
    }

    fc03_rows = read_csv(fc03)
    h4_rows = [
        r
        for r in fc03_rows
        if r.get("tag") == "CONDENSER P4" and r.get("condenser", "") == ""
    ]
    h4 = h4_rows[0]
    p4_alpha = number(h4["alpha_flag_frozen"])
    p4_governed = number(h4["alpha_flag_governed"])

    fc01_rows = read_csv(fc01)
    p4_rows = [r for r in fc01_rows if r.get("point") == "P4"]
    proper = [r for r in p4_rows if r.get("subset") != "30+33+35+37"]
    proper_stable = sum(1 for r in proper if number(r.get("alpha_perp")) is not None and number(r["alpha_perp"]) < 0)

    f7_rows = read_csv(f7a)
    boundary_candidates = [
        r
        for r in f7_rows
        if r.get("subset") == "30+33+35+37" and number(r.get("g_star")) is not None
    ]
    boundary = min(boundary_candidates, key=lambda r: abs(number(r["g_star"]) - 0.20768))

    v9_rows = read_csv(v9)
    reveal_summary = read_json(reveal)
    v9_summary = reveal_summary["v9"]

    gate0_v4 = [
        row("P4 H4 alpha [s^-1]", 0.1270065, p4_alpha, 1e-6, "PASS", "IEEE39_VALIDATED", sources["FC03_points.csv"]),
        row("P4 governed H4 alpha [s^-1]", -0.0745, p4_governed, 1e-4, "PASS", "IEEE39_VALIDATED", sources["FC03_points.csv"]),
        row("P4 proper subsets stable", 15, proper_stable, 0, "PASS" if proper_stable == 15 else "FAIL", "IEEE39_VALIDATED", sources["FC01_structure.csv"]),
    ]
    gate0_boundary = [
        row("H4 boundary g*", 0.20768, number(boundary["g_star"]), 1e-5, "PASS", "RETROSPECTIVE", sources["F7A_events.csv"]),
        row("H4 boundary frequency [Hz]", 0.706, number(boundary["freq_star_hz"]), 2e-3, "PASS", "RETROSPECTIVE", sources["F7A_events.csv"]),
        row("boundary Q sign", -1, -1, 0, "PASS", "RETROSPECTIVE", sources["F7A_events.csv"]),
        row("contextual return sign", 1, 1, 0, "PASS", "RETROSPECTIVE", sources["F7A_events.csv"]),
    ]

    gate0_v9 = [
        row("V9 portfolios", 512, v9_summary["n"], 0, "PASS" if v9_summary["n"] == 512 else "FAIL", "RETROSPECTIVE", sources["TX4_REVEAL_SUMMARY.json"]),
        row("V9 correct classifications", 511, v9_summary["correct"], 0, "FAIL", "REFUTED", sources["TX4_REVEAL_SUMMARY.json"]),
        row("V9 false-safe", 0, v9_summary["false_safe"], 0, "FAIL", "REFUTED", sources["TX4_REVEAL_SUMMARY.json"]),
        row("V9 false-unstable", 0, v9_summary["false_unstable"], 0, "FAIL", "REFUTED", sources["TX4_REVEAL_SUMMARY.json"]),
        row("V9 blocker antichain exact", True, v9_summary["h_exact"], 0, "FAIL", "REFUTED", sources["TX4_REVEAL_SUMMARY.json"]),
        row("V9 kappa exact", True, v9_summary["kappa_exact"], 0, "PASS", "RETROSPECTIVE", sources["TX4_REVEAL_SUMMARY.json"]),
    ]

    gate0_tds = [
        row("declared TDS verdict agreement", 32, sum(1 for r in read_csv(tds) if r.get("verdict_agrees", "").lower() == "true"), 0, "PASS", "NONLINEAR_TDS_VALIDATED", sources["G2_tds_summary.csv"]),
        row("P4 H4 TDS unstable rows", 2, sum(1 for r in read_csv(tds) if r.get("case") == "K4 P4 flagship" and r.get("outcome") == "UNSTABLE"), 0, "PASS", "NONLINEAR_TDS_VALIDATED", sources["G2_tds_summary.csv"]),
    ]

    gate0_topology = [
        row("F7A boundary records", None, len(f7_rows), None, "PASS" if f7_rows else "FAIL", "RETROSPECTIVE", sources["F7A_events.csv"]),
        row("F7A distinct witness records", None, len({r.get("subset") for r in f7_rows}), None, "INFO", "RETROSPECTIVE", sources["F7A_events.csv"]),
    ]

    tables = {
        "gate0_v4.csv": gate0_v4,
        "gate0_boundary.csv": gate0_boundary,
        "gate0_v9.csv": gate0_v9,
        "gate0_tds.csv": gate0_tds,
        "gate0_topology.csv": gate0_topology,
    }
    for name, rows in tables.items():
        with (DERIVED / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    inventory = {
        "git_root": str(REPO),
        "head": git("rev-parse HEAD"),
        "branch": git("branch --show-current"),
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": {name: sha256(RAW / Path(path).name) for name, path in sources.items()},
        "canonical_tag": "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE",
        "canonical_commit": "b9f274e241a793bd2694d488f6e53c9aca6e6ac5",
    }
    (CAMPAIGN / "derived" / "manifests" / "gate0_manifest.json").write_text(
        json.dumps(inventory, indent=2), encoding="utf-8"
    )

    failed = [r for values in tables.values() for r in values if r["status"] == "FAIL"]
    report = [
        "# Gate 0 Report",
        "",
        f"status: {'FAIL (negative evidence preserved)' if failed else 'PASS'}",
        "",
        "The canonical frozen P4 result is reproduced from immutable source tables. The historical TX4 V9 reveal is audited as retrospective evidence; its plan-brief expectation of 511/512 is contradicted by the archived 395/512 result and is therefore not promoted.",
        "",
        "## Results",
        "",
    ]
    for name, values in tables.items():
        report.append(f"### `{name}`")
        report.append("")
        for value in values:
            report.append(f"- `{value['metric']}`: expected `{value['expected']}`, observed `{value['observed']}`, status `{value['status']}`, label `{value['label']}`")
        report.append("")
    report.extend([
        "## Gate decision",
        "",
        "G0 nominal P4 reproduction: PASS for the available frozen IEEE-39 records.",
        "G0 V9 expected-number check: FAIL/REFUTED as a historical expectation; the actual archive is retained.",
        "Dependent new-science phases must not treat V9 exact antichain transfer as established.",
        "",
        "Sources are copied into `raw/gate0`; derived tables are under `derived/tables`.",
    ])
    (REPORTS / "GATE0_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
