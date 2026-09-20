"""Run final machine-checkable invariants for the corrected TX4 bundle."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "TX4_AUDIT_FIGURES"
PARAMS = ("g", "k", "t", "h", "epsilon", "damping", "inertia")
PORTFOLIOS = {
    "BASE", "30", "33", "30+33", "35", "30+35", "33+35", "30+33+35",
    "37", "30+37", "33+37", "30+33+37", "35+37", "30+35+37",
    "33+35+37", "30+33+35+37",
}


def row(checks, name, passed, observed, expected, note):
    checks.append(
        {
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
            "note": note,
        }
    )


def validate_campaign(checks: list[dict[str, object]], campaign: str, expected_conditions: int) -> pd.DataFrame:
    df = pd.read_parquet(RESULTS / f"TX4_{campaign}_ALL16_EXACT.parquet")
    row(checks, f"{campaign}_row_count", len(df) == expected_conditions * 16, len(df), expected_conditions * 16, "exact all-16 output row count")
    counts = df.groupby("condition_id").portfolio.nunique()
    row(checks, f"{campaign}_condition_count", len(counts) == expected_conditions, len(counts), expected_conditions, "unique exact conditions")
    row(checks, f"{campaign}_16_portfolios_each", bool((counts == 16).all()), int((counts == 16).sum()), expected_conditions, "every condition has all 16 portfolios")
    row(checks, f"{campaign}_portfolio_catalog", set(df.portfolio.unique()) == PORTFOLIOS, sorted(set(df.portfolio.unique())), sorted(PORTFOLIOS), "portfolio catalog exact")
    row(checks, f"{campaign}_exact_source", set(df.evaluation_source.dropna().unique()) == {"EXACT_DAE"}, sorted(df.evaluation_source.dropna().unique()), ["EXACT_DAE"], "fallback rows are exact DAE")
    truth = pd.read_csv(RESULTS / f"TX4_{campaign}_EXACT_BLOCKER_TRUTH.csv")
    row(checks, f"{campaign}_truth_condition_count", len(truth) == expected_conditions, len(truth), expected_conditions, "condition-level truth table")
    row(checks, f"{campaign}_implication_exact_present", bool((~truth.EXACT_H4 | truth.H4_PRESENT).all()), int((~truth.EXACT_H4 | truth.H4_PRESENT).sum()), expected_conditions, "EXACT_H4 implies H4_PRESENT")
    row(checks, f"{campaign}_implication_present_unstable", bool((~truth.H4_PRESENT | truth.H4_UNSTABLE).all()), int((~truth.H4_PRESENT | truth.H4_UNSTABLE).sum()), expected_conditions, "H4_PRESENT implies H4_UNSTABLE")
    row(checks, f"{campaign}_implication_exact_noncomposable", bool((~truth.EXACT_H4 | truth.NONCOMPOSABLE).all()), int((~truth.EXACT_H4 | truth.NONCOMPOSABLE).sum()), expected_conditions, "EXACT_H4 implies NONCOMPOSABLE")
    recomputed_delta = np.minimum(truth.h4_alpha_EM.to_numpy(float), -truth.max_proper_alpha_EM.to_numpy(float))
    row(checks, f"{campaign}_delta_formula", bool(np.allclose(recomputed_delta, truth.delta_H4_true_minimality.to_numpy(float), equal_nan=True)), "allclose", "allclose", "true minimality delta formula")
    for idx, truth_row in truth.iterrows():
        h0 = json.loads(truth_row.H0_json)
        sets = [frozenset() if x == "BASE" else frozenset(int(v) for v in x.split("+")) for x in h0]
        antichain = all(not (a < b or b < a) for i, a in enumerate(sets) for b in sets[i + 1 :])
        if not antichain:
            row(checks, f"{campaign}_h0_antichain", False, truth_row.condition_id, "all minimal", "strict subset found")
            break
    else:
        row(checks, f"{campaign}_h0_antichain", True, expected_conditions, expected_conditions, "every H0 is an antichain")
    return truth


def main() -> None:
    checks: list[dict[str, object]] = []
    q_truth = validate_campaign(checks, "QMC", 4096)
    m_truth = validate_campaign(checks, "MC", 5000)
    master = pd.read_parquet(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet")
    allowed = {"EXACT_DAE", "SURROGATE_EXTRATREES", "DERIVED_FROM_EXACT", "UNKNOWN"}
    row(checks, "master_evaluation_source_values", set(master.evaluation_source.unique()) <= allowed, sorted(master.evaluation_source.unique()), sorted(allowed), "registered provenance vocabulary")
    row(checks, "master_evaluation_source_row_count", len(master) == 329440, len(master), 329440, "every legacy master row carries provenance")
    for campaign, truth in (("QMC", q_truth), ("MC", m_truth)):
        exact = pd.read_parquet(RESULTS / f"TX4_{campaign}_ALL16_EXACT.parquet").drop_duplicates("condition_id")
        original = master[master.campaign == campaign].drop_duplicates("condition_id")
        merged = exact[["condition_id", *PARAMS]].merge(original[["condition_id", *PARAMS]], on="condition_id", suffixes=("_exact", "_master"))
        matched = len(merged) == len(exact) and all(np.allclose(merged[f"{p}_exact"], merged[f"{p}_master"], rtol=0, atol=0) for p in PARAMS)
        row(checks, f"{campaign}_coordinates_unchanged", matched, len(merged), len(exact), "exact fallback reuses retained coordinates bit-for-bit")
    eta = pd.read_csv(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv")
    row(checks, "eta_stratum_complete", len(eta) == 128 and (eta.status == "OK").all(), int((eta.status == "OK").sum()), 128, "128 exact stratified local/collective rows")
    row(checks, "eta_nonnegative", bool((eta.eta_corrected >= 0).all()), float(eta.eta_corrected.min()), ">=0", "singular-value eta is nonnegative")
    gate = json.loads((RESULTS / "TX4_SURROGATE_GATE.json").read_text(encoding="utf-8"))
    row(checks, "surrogate_gate_recorded_failure", gate["strict_gate_pass"] is False, gate["strict_gate_pass"], False, "exact fallback was required")
    row(checks, "exact_fallback_complete", (RESULTS / "TX4_AUDIT_EXACT_PROGRESS.json").exists() and json.loads((RESULTS / "TX4_AUDIT_EXACT_PROGRESS.json").read_text(encoding="utf-8"))["status"] == "COMPLETE", "COMPLETE", "COMPLETE", "checkpoint manifest")
    required = [
        RESULTS / "TX4_QMC_ALL16_EXACT.parquet", RESULTS / "TX4_MC_ALL16_EXACT.parquet",
        RESULTS / "TX4_AUDIT_H4_FLAG_COMPARISON.csv", RESULTS / "TX4_AUDIT_ETA_H4.csv",
        RESULTS / "TX4_ROBUSTNESS_DATA_PROVENANCE.csv", RESULTS / "TX4_SURROGATE_MINIMALITY_VALIDATION.csv",
        RESULTS / "TX4_SURROGATE_BLOCKER_ERRORS.csv", RESULTS / "TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv",
    ]
    row(checks, "required_audit_outputs", all(path.exists() for path in required), sum(path.exists() for path in required), len(required), "required audit files present")
    figures = sorted(FIGURES.glob("F*.png"))
    row(checks, "corrected_figure_count", len(figures) == 10, len(figures), 10, "F1-F10 corrected PNG set")
    all_pass = all(item["status"] == "PASS" for item in checks)
    pd.DataFrame(checks).to_csv(RESULTS / "TX4_ROBUSTNESS_INVARIANT_CHECKS.csv", index=False)
    print(json.dumps({"all_pass": all_pass, "n_checks": len(checks), "n_fail": sum(x["status"] == "FAIL" for x in checks)}, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

