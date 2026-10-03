#!/usr/bin/env python3
"""Audit poster-level arithmetic against the frozen, versioned source snapshot.

This checks counts and provenance used by the poster. It does not rerun the
power-flow, DAE, or eigenvalue simulations that produced the snapshot.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = HERE / "research/bnd_h4_mechanism/results/20260927T144629Z_d0fecb32_ias26_060_operating_v1"
F2 = HERE / "research/bnd_h4_mechanism/results/20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1/tables"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def yes(row: dict[str, str], key: str) -> bool:
    return row[key].strip().upper() == "TRUE"


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_manifest() -> None:
    manifest = json.loads((HERE / "EVIDENCE_MANIFEST.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        source = ROOT / entry["path"]
        # Git may check text files out with CRLF on Windows; hash LF-normalized
        # bytes so the frozen snapshot verifies on both Windows and Unix.
        digest = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        check(digest == entry["sha256"], f"Frozen source changed: {entry['path']}")


def check_macros(expected: dict[str, str]) -> None:
    source = (HERE / "generated/results.tex").read_text(encoding="utf-8")
    for name, value in expected.items():
        line = rf"\newcommand{{\{name}}}{{{value}}}"
        check(line in source, f"Generated macro does not match source: {name}")


def main() -> None:
    check_manifest()

    scenario = rows(RUN / "derived/SCENARIO_METRICS.csv")
    valid = [row for row in scenario if yes(row, "scenario_valid")]
    check(len(scenario) == 1000 and len(valid) == 967, "Scenario denominator changed")
    check(len({row["scenario_id"] for row in scenario}) == 1000, "Repeated scenario ID")
    no_blocker = sum(row["h4_class"] == "STABLE" for row in valid)
    h4_minimal = sum(yes(row, "h4_minimal_blocker_persists") for row in valid)
    triples = sum(row["h4_class"] == "UNSTABLE" and not yes(row, "h4_minimal_blocker_persists") for row in valid)
    collective = sum(yes(row, "phenomenon_persists") for row in valid)
    improved = sum(yes(row, "intervention_improves_h4") for row in valid)
    rescued = sum(yes(row, "intervention_rescues_h4") for row in valid)
    deteriorated = sum(yes(row, "intervention_deteriorates_h4") for row in valid)
    check((no_blocker, h4_minimal, triples, collective) == (412, 190, 365, 555), "Blocker classes changed")
    check(no_blocker + h4_minimal + triples == len(valid), "Blocker classes do not partition valid cases")
    check(h4_minimal + triples == collective, "Collective subtotal changed")
    check((improved, rescued, deteriorated) == (594, 186, 373), "Intervention outcomes changed")
    check(improved + deteriorated == len(valid), "Intervention outcomes do not partition valid cases")
    check(all(not yes(row, "intervention_rescues_h4") or yes(row, "intervention_improves_h4") for row in valid),
          "A rescue is not in the improvement group")

    event = {row["event"]: row for row in rows(RUN / "tables/MC_EVENT_RATES.csv")}
    for name, count in {
        "E1_H4_MINIMAL_BLOCKER_PERSISTS": h4_minimal,
        "E2_COLLECTIVE_PHENOMENON_PERSISTS": collective,
        "E3_INTERVENTION_IMPROVES": improved,
        "E4_INTERVENTION_RESCUES_H4": rescued,
        "E5_INTERVENTION_DETERIORATES": deteriorated,
    }.items():
        check(int(event[name]["numerator"]) == count and int(event[name]["denominator"]) == len(valid),
              f"Event table disagrees: {name}")

    census = rows(HERE / "research/results/FINAL_CLOSURE/FC10_census_lattice.csv")
    summary = json.loads((HERE / "research/results/FINAL_CLOSURE/FC10_summary.json").read_text(encoding="utf-8"))
    A = sum(yes(row, "A_final_stable") for row in census)
    B = sum(yes(row, "B_safe_path") for row in census)
    C = sum(yes(row, "C_any_order_safe") for row in census)
    check(len(census) == 512 and (A, B, C) == (327, 327, 324), "Path census changed")
    check((A, B, C) == (summary["A"], summary["B"], summary["C"]), "Path summary disagrees with census")
    check(sum(yes(row, "A_final_stable") and not yes(row, "C_any_order_safe") for row in census) == 3,
          "Order-dependent target count changed")

    baseline = {row["method_id"]: row for row in rows(HERE / "research/results/20260911_BASELINE_COMPARISON.csv")}
    for method, expected in {"B3": "22/52", "B4": "31/52", "B5": "38/52"}.items():
        check(baseline[method]["taskB_F10_52_H_exact"] == expected, f"Baseline {method} changed")

    policy = rows(HERE / "research/results/FINAL_CLOSURE/figures/FIG1_policy_map_source.csv")
    check(len(policy) == 120390 and len({row["H_perp"] for row in policy}) == 30,
          "Policy map A changed")
    claims = {row["id"]: row for row in rows(HERE / "research/results/20260911_CLAIM_MATRIX.csv")}
    check(re.search(r"\(30/36/16 distinct H", claims["B01"]["claim"]) is not None,
          "Three policy-sweep counts changed")

    boundary = rows(F2 / "F2_BOUNDARY_PORT_AUDIT.csv")
    check(len(boundary) == 4, "Boundary port count changed")
    check(max(float(row["g_root"]) for row in boundary) - min(float(row["g_root"]) for row in boundary) < 1e-6,
          "Boundary gain is inconsistent")
    check(abs(sum(float(row["g_root"]) for row in boundary) / 4 - 0.207681) < 1e-6,
          "Boundary gain changed")

    readiness = json.loads((RUN / "reports/IAS2026_POSTER_READINESS.json").read_text(encoding="utf-8"))
    tds = readiness["TDS"]
    check((tds["n_sign_agreement"], tds["n_completed_finite_primary_fits"], tds["n_planned"],
           tds["run_status_counts"]["SOLVER_FAILED"]) == (56, 60, 90, 30), "TDS counts changed")
    cross = rows(RUN / "tables/IAS26-FINAL_EVIDENCE_STRIP.csv")[0]
    check((cross["verdict_agreement"], cross["N_case_comparisons"]) == ("100", "100"),
          "Cross-code agreement changed")

    check_macros({
        "NumScenarios": "1000", "NumValidMC": "967", "NumInfeasibleMC": "33",
        "NumStableHFourScenarios": "412", "NoBlockerPercent": "42.6",
        "BlockerPersistenceCount": "190", "HFourMinimalPercent": "19.6",
        "NumAlternativeWitnessScenarios": "365", "AlternativeWitnessPercent": "37.7",
        "NumCollectiveBlockerScenarios": "555", "CollectiveBlockerPercent": "57.4",
        "InterventionImproves": "594", "InterventionRescues": "186",
        "InterventionDeteriorates": "373", "NumStableTargets": "327",
        "NumOrderDependentTargets": "3", "TDSAgreement": "56/60",
        "TDSSolverFailures": "30", "NumDistinctPolicyA": "30",
        "NumDistinctPolicyB": "36", "NumDistinctPolicyC": "16",
    })
    print("Poster evidence audit passed: frozen files, scenario partition, policy map, path census, baselines, boundary, and model checks.")


if __name__ == "__main__":
    main()
