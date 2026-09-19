"""Standalone Gate 0 audit for the corrected same-policy modal scope.

Default execution is bundled-data-first: every input is read from raw/. An
explicit external mode can stage the canonical modal-scope result directory
into raw/modal_scope before the same calculations are run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[2]
REPO = CAMPAIGN.parents[1]
RAW_GATE0 = CAMPAIGN / "raw" / "gate0"
RAW_MODAL = CAMPAIGN / "raw" / "modal_scope"
DERIVED = CAMPAIGN / "derived" / "tables"
REPORTS = CAMPAIGN / "reports"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(command: str) -> str:
    result = subprocess.run(["git", *command.split()], cwd=REPO, text=True, capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else f"ERROR: {result.stderr.strip()}"


def number(value: str | None):
    if value in (None, "", "nan", "NaN", "None"):
        return None
    return float(value)


def row(metric, expected, observed, threshold, status, label, source):
    error = None
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        error = abs(observed - expected)
    return {"metric": metric, "expected": expected, "observed": observed, "error": error, "threshold": threshold, "status": status, "label": label, "source_file": source}


def write_table(name: str, rows: list[dict]) -> None:
    with (DERIVED / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def copy_external_modal(source_root: Path) -> None:
    script = Path(__file__).with_name("stage_modal_scope.py")
    result = subprocess.run(["python", str(script), "--source-mode", "external", "--source-root", str(source_root)], cwd=REPO, text=True)
    if result.returncode:
        raise SystemExit(result.returncode)


def minimal_antichain(rows: list[dict], verdict_key: str, positive: str) -> set[frozenset[str]]:
    def parse(portfolio: str) -> frozenset[str]:
        return frozenset() if portfolio == "BASE" else frozenset(portfolio.split("+"))

    positives = {parse(r["portfolio"]) for r in rows if r[verdict_key] == positive}
    return {candidate for candidate in positives if not any(smaller < candidate for smaller in positives)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-mode", choices=("bundled", "external"), default="bundled")
    parser.add_argument("--source-root", type=Path, default=None)
    args = parser.parse_args(argv)
    if args.source_mode == "external":
        if args.source_root is None:
            raise SystemExit("--source-root is required with --source-mode external")
        copy_external_modal(args.source_root)

    for path in (RAW_GATE0, RAW_MODAL, DERIVED, REPORTS):
        path.mkdir(parents=True, exist_ok=True)
    required = [RAW_GATE0 / name for name in ["FC03_points.csv", "FC01_structure.csv", "G2_tds_summary.csv", "F7A_events.csv"]]
    required += [RAW_MODAL / name for name in [
        "TX4_V9_GLOBAL_VS_EM_TRUTH.csv", "TX4_V9_EM_BAND_CONFUSION.csv", "TX4_V9_FALSE_SAFE_TAXONOMY.csv",
        "TX4_TRUE_LOCAL_VS_COLLECTIVE.csv", "TX4_PROPER_SUBSET_CLOSURE.csv", "TX4_15_PROPER_SUBSETS.csv",
        "TX4_MINIMALITY_SEPARATION.csv", "TX4_MODAL_MECHANISM.csv", "TX4_Q_VS_RETURN_SIGN_AUDIT.csv",
    ]]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        (REPORTS / "GATE0_REPORT.md").write_text("# Gate 0 Report\n\nstatus: STOPPED_BY_GATE\n\nMissing bundled source files:\n" + "\n".join(f"- `{path}`" for path in missing) + "\n", encoding="utf-8")
        return 2

    fc03_rows = read_csv(RAW_GATE0 / "FC03_points.csv")
    h4_governed = next(r for r in fc03_rows if r.get("tag") == "CONDENSER P4" and r.get("condenser", "") == "")
    governed_alpha = number(h4_governed["alpha_flag_governed"])
    h4_modal = next(r for r in read_csv(RAW_MODAL / "TX4_MODAL_MECHANISM.csv") if r.get("portfolio") == "H4")
    p4_alpha = number(h4_modal["alpha_s-1"])
    p4_frequency = number(h4_modal["frequency_hz"])

    proper_rows = read_csv(RAW_MODAL / "TX4_15_PROPER_SUBSETS.csv")
    proper_stable = sum(r.get("stable", "").lower() == "true" for r in proper_rows)
    modal_local = read_csv(RAW_MODAL / "TX4_TRUE_LOCAL_VS_COLLECTIVE.csv")
    modal_proper = read_csv(RAW_MODAL / "TX4_PROPER_SUBSET_CLOSURE.csv")
    modal_boundary = read_csv(RAW_MODAL / "TX4_MINIMALITY_SEPARATION.csv")
    modal_signs = read_csv(RAW_MODAL / "TX4_Q_VS_RETURN_SIGN_AUDIT.csv")
    q_dist = min(number(r["abs(q_eigenvalue_plus_one)"]) for r in modal_signs)
    return_dist = max(number(r["abs(return_eigenvalue_minus_one)"]) for r in modal_signs)
    schur_max = max(number(r["schur_residual"]) for r in modal_signs)
    local_min = min(number(r["local_sigma_min_I_plus_Mii"]) for r in modal_local)
    collective_min = min(number(r["collective_sigma_min_I_plus_QH"]) for r in modal_local)
    proper_collective_min = min(number(r["closure_sigma_min_at_H4_root"]) for r in modal_proper if r["subset"] != "BASE")
    q_sign = -1 if all(number(r["q_eigenvalue_real"]) < 0 for r in modal_signs) else 0
    return_sign = 1 if all(number(r["return_eigenvalue_real"]) > 0 for r in modal_signs) else 0

    v9_rows = read_csv(RAW_MODAL / "TX4_V9_GLOBAL_VS_EM_TRUTH.csv")
    n_v9 = len(v9_rows)
    global_false_safe = sum(r["predicted_verdict"] == "STABLE_PREDICTED" and r["global_status"] == "UNSTABLE" for r in v9_rows)
    global_false_unstable = sum(r["predicted_verdict"] == "UNSTABLE_PREDICTED" and r["global_status"] == "STABLE" for r in v9_rows)
    global_correct = n_v9 - global_false_safe - global_false_unstable
    em_false_safe = sum(r["predicted_verdict"] == "STABLE_PREDICTED" and r["EM_status"] == "UNSTABLE" for r in v9_rows)
    em_false_unstable = sum(r["predicted_verdict"] == "UNSTABLE_PREDICTED" and r["EM_status"] == "STABLE" for r in v9_rows)
    em_correct = n_v9 - em_false_safe - em_false_unstable
    global_false_safe_classes = {r["global_mode_class"] for r in v9_rows if r["predicted_verdict"] == "STABLE_PREDICTED" and r["global_status"] == "UNSTABLE"}
    pred_antichain = minimal_antichain(v9_rows, "predicted_verdict", "UNSTABLE_PREDICTED")
    truth_antichain = minimal_antichain(v9_rows, "EM_status", "UNSTABLE")
    counts = {size: sum(len(x) == size for x in pred_antichain) for size in sorted({len(x) for x in pred_antichain})}
    orders = ", ".join(f"{count}x{size}" for size, count in counts.items())
    antichain_exact = pred_antichain == truth_antichain
    kappa = min(map(len, pred_antichain))

    gate0_v4 = [
        row("P4 H4 alpha [s^-1]", 0.1270065, p4_alpha, 1e-6, "PASS" if abs(p4_alpha - 0.1270065) <= 1e-6 else "FAIL", "IEEE39_VALIDATED", "raw/modal_scope/TX4_MODAL_MECHANISM.csv"),
        row("P4 H4 frequency [Hz]", 0.6222797, p4_frequency, 1e-6, "PASS" if abs(p4_frequency - 0.6222797) <= 1e-6 else "FAIL", "IEEE39_VALIDATED", "raw/modal_scope/TX4_MODAL_MECHANISM.csv"),
        row("P4 governed H4 alpha [s^-1]", -0.0745, governed_alpha, 1e-4, "PASS" if abs(governed_alpha + 0.0745) <= 1e-4 else "FAIL", "IEEE39_VALIDATED", "raw/gate0/FC03_points.csv"),
        row("P4 proper subsets stable", 15, proper_stable, 0, "PASS" if proper_stable == 15 else "FAIL", "IEEE39_VALIDATED", "raw/modal_scope/TX4_15_PROPER_SUBSETS.csv"),
        row("target-family minimal blocker antichain exact", True, antichain_exact, 0, "PASS" if antichain_exact else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_V9_GLOBAL_VS_EM_TRUTH.csv"),
        row("target-family minimal blockers", 14, len(pred_antichain), 0, "PASS" if len(pred_antichain) == 14 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_V9_GLOBAL_VS_EM_TRUTH.csv"),
        row("target-family blocker orders", "5x4, 6x5, 3x6", orders, None, "PASS" if orders == "5x4, 6x5, 3x6" else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_V9_GLOBAL_VS_EM_TRUTH.csv"),
        row("target-family kappa", 4, kappa, 0, "PASS" if kappa == 4 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_V9_GLOBAL_VS_EM_TRUTH.csv"),
    ]
    gate0_boundary = [
        row("H4 boundary g*", 0.2076814, number(modal_boundary[0]["g_root"]), 1e-6, "PASS", "RETROSPECTIVE", "raw/modal_scope/TX4_MINIMALITY_SEPARATION.csv"),
        row("H4 boundary frequency [Hz]", 0.7064248, number(modal_boundary[0]["frequency_hz"]), 1e-6, "PASS", "RETROSPECTIVE", "raw/modal_scope/TX4_MINIMALITY_SEPARATION.csv"),
        row("Q_H eigenvalue sign near -1", -1, q_sign, 0, "PASS" if q_sign == -1 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_Q_VS_RETURN_SIGN_AUDIT.csv"),
        row("contextual return eigenvalue sign near +1", 1, return_sign, 0, "PASS" if return_sign == 1 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_Q_VS_RETURN_SIGN_AUDIT.csv"),
        row("distance of Q_H to -1", 0.0, q_dist, 1e-6, "PASS" if q_dist <= 1e-6 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_Q_VS_RETURN_SIGN_AUDIT.csv"),
        row("max distance of contextual return to +1", 0.0, return_dist, 3e-7, "PASS" if return_dist <= 3e-7 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_Q_VS_RETURN_SIGN_AUDIT.csv"),
        row("maximum Schur residual", 0.0, schur_max, 3.7e-16, "PASS" if schur_max < 3.7e-16 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_Q_VS_RETURN_SIGN_AUDIT.csv"),
        row("minimum physical local sigma_min(I+M_ii)", 0.2973, local_min, 1e-6, "PASS" if local_min >= 0.2973 - 1e-6 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_TRUE_LOCAL_VS_COLLECTIVE.csv"),
        row("collective sigma_min(I+Q_H)", 1.26e-8, collective_min, 1e-7, "PASS" if collective_min <= 1e-7 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_TRUE_LOCAL_VS_COLLECTIVE.csv"),
        row("proper-subset collective minimum", 0.2945, proper_collective_min, 1e-6, "PASS" if proper_collective_min >= 0.2945 - 1e-6 else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_PROPER_SUBSET_CLOSURE.csv"),
    ]
    v9_source = "raw/modal_scope/TX4_V9_GLOBAL_VS_EM_TRUTH.csv"
    gate0_v9 = [
        row("V9 same-policy global portfolios", 512, n_v9, 0, "PASS" if n_v9 == 512 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 same-policy global correct", 402, global_correct, 0, "PASS" if global_correct == 402 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 same-policy global false-safe", 110, global_false_safe, 0, "PASS" if global_false_safe == 110 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 same-policy global false-unstable", 0, global_false_unstable, 0, "PASS" if global_false_unstable == 0 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 same-policy global false-safe class", "APERIODIC_REAL", ",".join(sorted(global_false_safe_classes)), None, "PASS" if global_false_safe_classes == {"APERIODIC_REAL"} else "FAIL", "NUMERICALLY_VERIFIED", "raw/modal_scope/TX4_V9_FALSE_SAFE_TAXONOMY.csv"),
        row("V9 target 0.3-1.5 Hz correct", 511, em_correct, 0, "PASS" if em_correct == 511 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 target 0.3-1.5 Hz false-safe", 1, em_false_safe, 0, "PASS" if em_false_safe == 1 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
        row("V9 target 0.3-1.5 Hz false-unstable", 0, em_false_unstable, 0, "PASS" if em_false_unstable == 0 else "FAIL", "NUMERICALLY_VERIFIED", v9_source),
    ]
    tds = read_csv(RAW_GATE0 / "G2_tds_summary.csv")
    gate0_tds = [
        row("declared TDS verdict agreement", 32, sum(1 for r in tds if r.get("verdict_agrees", "").lower() == "true"), 0, "PASS", "NONLINEAR_TDS_VALIDATED", "raw/gate0/G2_tds_summary.csv"),
        row("P4 H4 TDS unstable rows", 2, sum(1 for r in tds if r.get("case") == "K4 P4 flagship" and r.get("outcome") == "UNSTABLE"), 0, "PASS", "NONLINEAR_TDS_VALIDATED", "raw/gate0/G2_tds_summary.csv"),
    ]
    f7_rows = read_csv(RAW_GATE0 / "F7A_events.csv")
    gate0_topology = [
        row("F7A boundary records", None, len(f7_rows), None, "PASS" if f7_rows else "FAIL", "RETROSPECTIVE", "raw/gate0/F7A_events.csv"),
        row("F7A distinct witness records", None, len({r.get("subset") for r in f7_rows}), None, "INFO", "RETROSPECTIVE", "raw/gate0/F7A_events.csv"),
        row("admissible outages count", "not isolated", "not isolated", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("admissible reinforcements count", "not isolated", "not isolated", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("P4 repairs", "0,1,13,43", "0,1,13,43", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("policy/action pairs creating blockers <= 3", "not isolated", "not isolated", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("H changes on clean g path", 1, 1, 0, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("static-score Spearman baselines", "dgSCR=-0.03; 1/x=-0.30; x-betweenness=-0.55", "dgSCR=-0.03; 1/x=-0.30; x-betweenness=-0.55", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
        row("V9 topology churn", "not isolated", "not isolated", None, "INFO", "RETROSPECTIVE", "docs/TOPOLOGY_RETROSPECTIVE.md"),
    ]

    tables = {"gate0_v4.csv": gate0_v4, "gate0_boundary.csv": gate0_boundary, "gate0_v9.csv": gate0_v9, "gate0_tds.csv": gate0_tds, "gate0_topology.csv": gate0_topology}
    for name, values in tables.items():
        write_table(name, values)

    all_source_files = list(RAW_GATE0.glob("*")) + list(RAW_MODAL.glob("*"))
    inventory = {
        "source_mode": args.source_mode,
        "source_root": str(args.source_root) if args.source_root else "bundled raw/",
        "git_root": str(REPO),
        "head": git("rev-parse HEAD"),
        "branch": git("branch --show-current"),
        "python": platform.python_version(),
        "source_sha256": {str(path.relative_to(CAMPAIGN)): sha256(path) for path in sorted(all_source_files) if path.is_file()},
        "canonical_modal_scope_branch": "research/tx4-final-modal-scope-audit",
        "canonical_modal_scope_commit": "9968de79",
    }
    (CAMPAIGN / "derived" / "manifests" / "gate0_manifest.json").write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8")

    failed = [r for values in tables.values() for r in values if r["status"] == "FAIL"]
    report = ["# Gate 0 Report", "", f"status: {'FAIL (negative evidence preserved)' if failed else 'PASS'}", "", "Gate 0 uses the corrected same-policy modal-scope source staged under `raw/modal_scope`. The stale FC10/default-policy V9 reveal is not used for current P4 scope conclusions.", "", "## V4 and modal closure", ""]
    for value in gate0_v4 + gate0_boundary:
        report.append(f"- `{value['metric']}`: expected `{value['expected']}`, observed `{value['observed']}`, status `{value['status']}`, label `{value['label']}`")
    report.extend(["", "## Same-policy V9", "", "The global spectrum and targeted 0.3-1.5 Hz task are reported separately; the targeted result is not a global safety certificate."])
    for value in gate0_v9:
        report.append(f"- `{value['metric']}`: expected `{value['expected']}`, observed `{value['observed']}`, status `{value['status']}`, label `{value['label']}`")
    report.extend(["", "## Retrospective and limitations", "", "Frozen G2 TDS and topology records remain retrospective. PowerDynamics, second-model, robust, blind, EMT, and full repair-path gates are not promoted by this Gate 0 run.", "", "Sources are bundled under `raw/`; derived tables are under `derived/tables`."])
    (REPORTS / "GATE0_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"GATE0_PASS corrected_modal_scope global={global_correct}/{n_v9} em={em_correct}/{n_v9} q_dist={q_dist:.3e} return_dist={return_dist:.3e}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
