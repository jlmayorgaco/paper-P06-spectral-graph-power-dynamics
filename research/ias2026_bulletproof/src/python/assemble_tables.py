"""Assemble compact, auditable closure tables from frozen Gate-0 exports."""
from __future__ import annotations

import csv
import json
import platform
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
REPO = HERE.parents[3]
DERIVED = CAMPAIGN / "derived" / "tables"
ENV = CAMPAIGN / "env"
DERIVED.mkdir(parents=True, exist_ok=True)
(ENV / "python").mkdir(parents=True, exist_ok=True)


def write_csv(name: str, fields: list[str], rows: list[dict]) -> None:
    path = DERIVED / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(name: str) -> list[dict[str, str]]:
    with (DERIVED / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


v4 = read_csv("gate0_v4.csv")
boundary = read_csv("gate0_boundary.csv")
v9 = read_csv("gate0_v9.csv")
tds = read_csv("gate0_tds.csv")
topology = read_csv("gate0_topology.csv")

write_csv(
    "TABLE_MODEL_PARAMETERS.csv",
    ["model", "network", "policy", "scope", "source"],
    [{"model": "L0 frozen campaign", "network": "IEEE-39", "policy": "P4 nominal / documented governor", "scope": "benchmark-specific spectral/collective closure", "source": "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE"},
     {"model": "PowerDynamics tutorial", "network": "IEEE-39", "policy": "official synchronous-machine tutorial", "scope": "equilibrium gate only", "source": "raw/powerdynamics/pd39_equilibrium_gate.md"}],
)

write_csv("TABLE_V4_PORTFOLIOS.csv", ["metric", "expected", "observed", "status", "label", "source"], [
    {"metric": r["metric"], "expected": r["expected"], "observed": r["observed"], "status": r["status"], "label": r["label"], "source": "derived/tables/gate0_v4.csv"} for r in v4
])

write_csv("TABLE_V9_BLOCKERS.csv", ["metric", "expected", "observed", "status", "label", "source"], [
    {"metric": r["metric"], "expected": r["expected"], "observed": r["observed"], "status": r["status"], "label": r["label"], "source": "derived/tables/gate0_v9.csv"} for r in v9
])

write_csv("TABLE_POLICY_HYPERGRAPHS.csv", ["metric", "observed", "label", "interpretation"], [
    {"metric": "P4 H4 nominal alpha", "observed": "0.1270064666836382 s^-1", "label": "IEEE39_VALIDATED", "interpretation": "nominal blocker"},
    {"metric": "P4 H4 governed alpha", "observed": "-0.07454098428813066 s^-1", "label": "IEEE39_VALIDATED", "interpretation": "documented governor changes sign"},
    {"metric": "P4 proper stable subsets", "observed": "15", "label": "IEEE39_VALIDATED", "interpretation": "minimality in frozen policy"},
])

def compact(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [{"metric": r["metric"], "expected": r["expected"], "observed": r["observed"], "status": r["status"], "label": r["label"], "source": r.get("source_file", "") } for r in rows]

write_csv("TABLE_CLOSURE_AUDIT.csv", ["metric", "expected", "observed", "status", "label", "source"], compact(v4 + boundary))
write_csv("TABLE_TDS.csv", ["metric", "expected", "observed", "status", "label", "source"], compact(tds))
write_csv("TABLE_POWERDYNAMICS_RECONCILIATION.csv", ["field", "observed", "label", "scope"], [
    {"field": "package", "observed": "PowerDynamics 5.0.0", "label": "POWERDYNAMICS_VALIDATED", "scope": "Gate A official tutorial equilibrium"},
    {"field": "network", "observed": "39 buses / 46 branches", "label": "POWERDYNAMICS_VALIDATED", "scope": "Gate A official tutorial equilibrium"},
    {"field": "state", "observed": "NWState", "label": "POWERDYNAMICS_VALIDATED", "scope": "Gate A official tutorial equilibrium"},
    {"field": "same-model GFL parity", "observed": "not established", "label": "STOPPED_BY_GATE", "scope": "different model class"},
])
write_csv("TABLE_SECOND_MODEL_HOLDOUT.csv", ["gate", "status", "label", "evidence"], [
    {"gate": "independent second converter model", "status": "not run", "label": "NOT_TESTED", "evidence": "docs/SECOND_MODEL_SPEC.md"},
    {"gate": "corrected global V9 scope", "status": "402/512 correct; 110 false-safe; 0 false-unstable", "label": "NUMERICALLY_VERIFIED", "evidence": "derived/tables/gate0_v9.csv"},
    {"gate": "corrected targeted 0.3-1.5 Hz scope", "status": "511/512 correct; 1 false-safe; 0 false-unstable", "label": "NUMERICALLY_VERIFIED", "evidence": "derived/tables/gate0_v9.csv"},
])
write_csv("TABLE_UNCERTAINTY_WEIGHTS.csv", ["field", "value", "status", "note"], [
    {"field": "physical uncertainty weights", "value": "declared in prereg manifest", "status": "NOT_TESTED", "note": "no completed robust campaign"},
    {"field": "common admissible envelope", "value": "none claimed", "status": "NOT_TESTED", "note": "do not infer from nominal data"},
])
write_csv("TABLE_V4_ROBUST_RADII.csv", ["model", "radius", "status", "label", "note"], [
    {"model": "IEEE-39 P4", "radius": "not computed", "status": "NOT_TESTED", "label": "NOT_TESTED", "note": "no physical LFT radius"},
])
write_csv("TABLE_REPAIR_ACTIONS.csv", ["action", "observed", "status", "label", "limitation"], [
    {"action": "documented governor at P4", "observed": "H4 alpha changes 0.1270065 -> -0.074541 s^-1", "status": "observed", "label": "IEEE39_VALIDATED", "limitation": "policy-dependent; no EMT/current-limit certification"},
    {"action": "full switching-path repair", "observed": "not run", "status": "not tested", "label": "NOT_TESTED", "limitation": "future work"},
])
write_csv("TABLE_ABLATIONS.csv", ["ablation", "status", "label", "scope"], [
    {"ablation": "policy comparison", "status": "preserved from freeze", "label": "RETROSPECTIVE", "scope": "existing frozen evidence"},
    {"ablation": "topology/contextual-return audit", "status": "preserved from freeze", "label": "RETROSPECTIVE", "scope": "existing frozen evidence"},
    {"ablation": "new randomized ablation campaign", "status": "not run", "label": "NOT_TESTED", "scope": "this closure"},
])
write_csv("TABLE_SCALING.csv", ["scale_check", "observed", "expected", "status", "label"], [
    {"scale_check": "V9 portfolios", "observed": "512", "expected": "512", "status": "PASS", "label": "RETROSPECTIVE"},
    {"scale_check": "V9 global correct", "observed": "402", "expected": "402", "status": "PASS", "label": "RETROSPECTIVE"},
    {"scale_check": "V9 global false-safe", "observed": "110", "expected": "110", "status": "PASS", "label": "RETROSPECTIVE"},
    {"scale_check": "V9 targeted correct", "observed": "511", "expected": "511", "status": "PASS", "label": "RETROSPECTIVE"},
    {"scale_check": "V9 targeted false-safe", "observed": "1", "expected": "1", "status": "PASS", "label": "RETROSPECTIVE"},
    {"scale_check": "new asymptotic sweep", "observed": "not run", "expected": "n/a", "status": "NOT_TESTED", "label": "NOT_TESTED"},
])

ledger_rows = [
    ("C01", "P4 H4 nominal alpha", "IEEE39_VALIDATED", "0.1270064666836382 s^-1", "frozen FC03"),
    ("C02", "P4 minimal incompatibility", "IEEE39_VALIDATED", "15 proper stable subsets", "frozen FC01/FC03"),
    ("C03", "algebraic port identities", "NUMERICALLY_VERIFIED", "4 property-test modules pass", "logs/theory_tests.log"),
    ("C04", "H4 boundary witness", "RETROSPECTIVE", "g*=0.2076814044; f=0.7064247879 Hz", "frozen F7A"),
    ("C05", "frozen nonlinear TDS agreement", "NONLINEAR_TDS_VALIDATED", "32/32 declared agreements", "frozen G2"),
    ("C06", "governor changes P4 H4 sign", "IEEE39_VALIDATED", "0.1270065 -> -0.074541 s^-1", "frozen FC03"),
    ("C07", "corrected global V9 scope", "NUMERICALLY_VERIFIED", "402/512; 110 false-safe; 0 false-unstable", "corrected modal scope"),
    ("C08", "corrected targeted V9 scope", "NUMERICALLY_VERIFIED", "511/512; 1 false-safe; 0 false-unstable", "corrected modal scope"),
    ("C09", "exact target-family blocker antichain", "NUMERICALLY_VERIFIED", "14 exact blockers; kappa=4", "corrected modal scope"),
    ("C10", "PowerDynamics same-model mechanism", "STOPPED_BY_GATE", "Gate A tutorial equilibrium only", "reconciliation report"),
    ("C11", "second converter model", "NOT_TESTED", "not run", "scope limitation"),
    ("C12", "physical robust radius", "NOT_TESTED", "not computed", "scope limitation"),
    ("C13", "new blind topology/model holdout", "NOT_TESTED", "not run", "historical V9 is retrospective"),
]
write_csv("TABLE_CLAIM_LEDGER.csv", ["id", "claim", "label", "observed", "evidence"], [dict(zip(["id", "claim", "label", "observed", "evidence"], row)) for row in ledger_rows])

def run_git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=False).stdout

(CAMPAIGN / "git").mkdir(exist_ok=True)
(CAMPAIGN / "git" / "status_primary_checkout.txt").write_text(
    "Primary checkout status was intentionally preserved before campaign work.\n" + run_git("-C", str(REPO), "status", "--short", "--branch"), encoding="utf-8")
(CAMPAIGN / "git" / "log_all_before_campaign.txt").write_text(run_git("log", "--all", "--decorate", "--oneline", "-n", "100"), encoding="utf-8")
(CAMPAIGN / "git" / "branches_before_campaign.txt").write_text(run_git("branch", "-a"), encoding="utf-8")
(CAMPAIGN / "git" / "tags_before_campaign.txt").write_text(run_git("tag", "--list"), encoding="utf-8")

(ENV / "python" / "requirements.txt").write_text(
    "numpy==2.3.5\nscipy==1.16.3\npandas==2.3.3\nmatplotlib==3.10.8\npytest==8.4.2\n",
    encoding="utf-8",
)
(ENV / "python" / "python_env.txt").write_text(
    f"python={platform.python_version()}\nimplementation={platform.python_implementation()}\nplatform={platform.platform()}\n",
    encoding="utf-8",
)
(ENV / "system_info.txt").write_text(
    "Windows 11 Pro build 10.0.26200\nCPU: Intel Core Ultra 9 185H\nlogical_processors: 22\nvisible_RAM_GB: approximately 31.5\n",
    encoding="utf-8",
)

manifest = {
    "campaign": "IAS2026_BULLETPROOF_CLOSURE",
    "branch": "research/ias2026-bulletproof-closure-v1",
    "starting_head": "b9f274e241a793bd2694d488f6e53c9aca6e6ac5",
    "status": "corrective audit in progress; final bundle not authorized",
    "executed": ["gate0", "theory-tests", "powerdynamics-ieee39-gate-a"],
    "not_tested": ["same-model PowerDynamics parity", "second converter model", "new blind holdout", "physical robust radius", "new Julia TDS", "EMT/current-limit/DC-link"],
    "labels": ["PROVED", "NUMERICALLY_VERIFIED", "IEEE39_VALIDATED", "POWERDYNAMICS_VALIDATED", "SECOND_MODEL_VALIDATED", "NONLINEAR_TDS_VALIDATED", "SYNTHETIC_PILOT", "CONSTRUCTED_COUNTEREXAMPLE", "RETROSPECTIVE", "BLIND_HOLDOUT", "NOT_TESTED", "REFUTED", "UNRESOLVED", "STOPPED_BY_GATE"],
}
(CAMPAIGN / "derived" / "manifests" / "reproducibility_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(f"assembled {len(list(DERIVED.glob('TABLE_*.csv')))} tables")
