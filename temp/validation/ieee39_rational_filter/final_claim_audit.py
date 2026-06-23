"""Create a final claim-readiness audit for the paper package.

This script intentionally separates mathematical/reduced-model evidence from
full ANDES IBR benchmark evidence. It does not rerun simulations; it audits the
versioned outputs produced by the validation scripts and writes a reviewer-facing
status report.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "final_claim_audit"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def pct(x: float) -> str:
    return f"{100.0 * x:.3g}%"


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def csv_bool(value: str) -> bool:
    return str(value).strip().lower() == "true"


def main() -> None:
    surrogate_path = ROOT / "outputs" / "ieee39_rational_filter" / "surrogate_summary.json"
    ibr_path = ROOT / "outputs" / "ieee39_ibr_cases" / "case_generation_summary.json"
    diag_path = ROOT / "outputs" / "ieee39_mode_diagnostics" / "baseline_positive_mode_diagnostics.json"
    phase1_sweep_path = ROOT / "outputs" / "ieee39_phase1_calibration" / "phase1_calibration_sweep.csv"
    phase1_best_path = ROOT / "outputs" / "ieee39_phase1_calibration" / "phase1_best_by_case.csv"
    physical_bridge_path = ROOT / "outputs" / "phase0_physical_bridge" / "phase0_physical_bridge_status.json"

    surrogate = load_json(surrogate_path)
    ibr_cases = load_json(ibr_path)
    diagnostics = load_json(diag_path)
    phase1_sweep = load_csv(phase1_sweep_path) if phase1_sweep_path.exists() else []
    phase1_best = load_csv(phase1_best_path) if phase1_best_path.exists() else []
    physical_bridge = load_json(physical_bridge_path) if physical_bridge_path.exists() else {}

    errors = surrogate["errors"]
    diag_median = errors["diag"]["median"]
    second_median = errors["second"]["median"]
    rqep_p95 = errors["rqep"]["p95"]
    second_p95 = errors["second"]["p95"]
    adaptive_p95 = errors["adaptive"]["p95"]

    n_ibr_cases = len(ibr_cases)
    n_ibr_ok = sum(1 for item in ibr_cases if item["test"]["ok"])
    n_ibr_stable = sum(1 for item in ibr_cases if item["test"]["poles"]["small_signal_stable"])
    has_gfl = any(item["test"]["inventory"]["REGCP1"] > 0 and item["test"]["inventory"]["PLL1"] > 0 for item in ibr_cases)
    has_gfm = any(item["test"]["inventory"]["REGF1"] > 0 for item in ibr_cases)
    n_phase1_runs = len(phase1_sweep)
    n_phase1_stable = sum(1 for row in phase1_sweep if csv_bool(row.get("small_signal_stable", "")))
    n_phase1_controller_stable = sum(1 for row in phase1_best if csv_bool(row.get("controller_stable", "")))
    n_phase1_diagnostic_stable = sum(1 for row in phase1_best if csv_bool(row.get("diagnostic_stable", "")))
    physical_status = physical_bridge.get("status", "not run")
    physical_results = physical_bridge.get("results", [])
    physical_base = physical_results[0] if physical_results else {}
    physical_crit = physical_base.get("critical_match", {}) if physical_base else {}
    physical_freq_error = physical_crit.get("relative_freq_error")
    physical_zeta_error = physical_crit.get("relative_zeta_error")

    base_diag = next(item for item in diagnostics if item["variant"] == "base")
    no_pss_diag = next(item for item in diagnostics if item["variant"] == "no_pss")

    braess_count = surrogate["braess_like_count"]
    n_lines = surrogate["n_lines_tested"]
    inertia_counts_path = ROOT / "outputs" / "inertia_reversal" / "inertia_reversal_counts.json"
    inertia_counts = load_json(inertia_counts_path) if inertia_counts_path.exists() else {}
    fixed_ratio = inertia_counts.get("fixed_damping_per_inertia", {})

    rows = [
        {
            "claim": "C1 damping-region coordinate system",
            "what_was_verified": "Exact in commuting/proportional damping by simultaneous diagonalization; figure shows non-proportional screen failure mode.",
            "control_or_comparison": "Commuting case versus non-proportional coupled-QEP case.",
            "not_proven": "The region language itself is not proven novel; outside commuting damping it is a screen, not a certificate.",
            "status": "theoretical/diagnostic",
            "confidence": "high for the mathematics; medium for novelty and planning value",
        },
        {
            "claim": "C2 closed-form second-order damping-margin correction",
            "what_was_verified": (
                "On the IEEE39-derived surrogate, median relative error improves "
                f"from {pct(diag_median)} diagonal to {pct(second_median)} second order; "
                "reduced PLL Monte Carlo reported in the manuscript improves 2.16% to 0.09% median error."
            ),
            "control_or_comparison": "Full QEP of the surrogate; diagonal all-mode screen; reduced QEP fallback.",
            "not_proven": "Not yet validated against full calibrated ANDES IBR poles; second order does not dominate large reduced QEP in tails.",
            "status": "verified on reduced models; benchmark pending",
            "confidence": "high for reduced models; medium for full IBR generalization",
        },
        {
            "claim": "C2 adaptive workflow diagonal -> second order -> reduced QEP",
            "what_was_verified": (
                "Tail behavior requires fallback: second-order p95 is "
                f"{pct(second_p95)}, reduced-QEP p95 is {pct(rqep_p95)}, "
                f"adaptive p95 is {pct(adaptive_p95)} on the IEEE39-derived surrogate."
            ),
            "control_or_comparison": "Always-diagonal, always-second-order, and reduced-QEP estimates.",
            "not_proven": "No universal trigger threshold is proved; false-safe behavior must be measured on calibrated ANDES IBR cases.",
            "status": "diagnostic workflow",
            "confidence": "medium",
        },
        {
            "claim": "C3 SG-to-IBR conversion ranking",
            "what_was_verified": (
                f"{n_ibr_ok}/{n_ibr_cases} generated SG-to-IBR ANDES cases load, solve power flow, and run eigenanalysis; "
                f"{n_ibr_stable}/{n_ibr_cases} is small-signal stable with generic uncalibrated gains. "
                f"Phase-1 ran {n_phase1_runs} case/profile audits; {n_phase1_controller_stable}/{n_ibr_cases} cases have a stable controller-only profile and "
                f"{n_phase1_diagnostic_stable}/{n_ibr_cases} have a stable diagnostic retained-SG-control-removal profile."
            ),
            "control_or_comparison": "ANDES pflow/eigenanalysis feasibility, device inventory checks, controller-only profiles, and diagnostic retained-control removal.",
            "not_proven": "No low-regret conversion ranking against full eigensolve has been demonstrated yet; stable cases are still calibration candidates, not final benchmark tuning.",
            "status": "case-construction and calibration audit",
            "confidence": "high for construction; low for planning superiority",
        },
        {
            "claim": "C4 damping-aware weak links and Braess-like reversals",
            "what_was_verified": (
                f"On the IEEE39-derived surrogate, {braess_count}/{n_lines} line reinforcements increase critical modal stiffness "
                "while reducing damping margin; the best damping-aware reinforcement is line 25-37."
            ),
            "control_or_comparison": "Frequency-only line effect versus finite post-action damping-margin change.",
            "not_proven": "Not yet a full ANDES IBR topology claim; generic QEP sensitivity prior art must be distinguished.",
            "status": "surrogate planning diagnostic",
            "confidence": "medium",
        },
        {
            "claim": "C5 inertia-placement sign reversal",
            "what_was_verified": (
                "With fixed local damping-per-inertia, "
                f"{fixed_ratio.get('both', 'NA')}/{fixed_ratio.get('tested', 'NA')} perturbations reduced both the "
                "diagonal modal screen and the full-QEP damping margin; the strongest reported example matches the analytic sign against finite differences."
            ),
            "control_or_comparison": "Fixed D/M nontrivial check versus fixed-D direct damping-dilution control.",
            "not_proven": "Not validated on calibrated ANDES IBR dynamics; the theorem is a reduced-model modal sign condition, not a universal inertia-planning law.",
            "status": "verified reduced-model sign diagnostic",
            "confidence": "medium",
        },
        {
            "claim": "C6 joint inertia-damping cross term",
            "what_was_verified": "The manuscript derives the local off-diagonal coupling term for simultaneous inertia and damping perturbations and shows that it reduces to the existing second-order damping correction when Delta M is zero.",
            "control_or_comparison": "Special-case reduction to the fixed-inertia second-order correction; separation of first-order diagonal shift from second-order modal coupling.",
            "not_proven": "No Monte Carlo or ANDES validation is claimed for finite SG-to-IBR conversion accuracy; large conversions must be recomputed.",
            "status": "theoretical local mechanism",
            "confidence": "medium for derivation; low for finite-conversion prediction",
        },
        {
            "claim": "ANDES IEEE39 benchmark readiness",
            "what_was_verified": (
                f"Packaged baseline has {base_diag['n_positive']} positive-real non-oscillatory modes; removing IEEEST gives "
                f"{no_pss_diag['n_positive']} positive modes while keeping the 1.37 Hz oscillatory pair essentially unchanged. "
                f"GFL/PLL cases exist: {yes_no(has_gfl)}; GFM cases exist: {yes_no(has_gfm)}. "
                f"Phase-1 found {n_phase1_stable}/{n_phase1_runs} stable case/profile runs. "
                f"The strict physical bridge audit status is {physical_status}"
                + (
                    f" (base critical frequency error {pct(float(physical_freq_error))}, damping-ratio error {pct(float(physical_zeta_error))})."
                    if physical_freq_error is not None and physical_zeta_error is not None
                    else "."
                )
            ),
            "control_or_comparison": "Baseline versus removed-model-family diagnostics; IBR device inventory checks.",
            "not_proven": "A defensible calibrated dynamic library and physically matched extraction of reduced L,M,D/control coupling from full ANDES DAE are still pending; the strict physical bridge audit remains blocked.",
            "status": "benchmark harness and calibration audit ready; final validation pending",
            "confidence": "high for harness status; low for final empirical claim",
        },
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "claim_audit.csv"
    json_path = OUT / "claim_audit.json"
    md_path = OUT / "final_readiness_report.md"

    fieldnames = ["claim", "what_was_verified", "control_or_comparison", "not_proven", "status", "confidence"]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    lines = [
        "# Final Claim-Readiness Audit",
        "",
        "This report is generated from versioned validation outputs. It does not rerun simulations.",
        "",
        "## Bottom Line",
        "",
        "The paper is ready as a rigorous methodological/master-thesis manuscript package. It is not yet ready to claim final IEEE Transactions empirical validation on full calibrated IBR dynamics.",
        "",
        "## Claim Status",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"### {row['claim']}",
                "",
                f"- Verified: {row['what_was_verified']}",
                f"- Control/comparison: {row['control_or_comparison']}",
                f"- Not proven: {row['not_proven']}",
                f"- Status: {row['status']}",
                f"- Confidence: {row['confidence']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Submission Decision",
            "",
            "Use the current manuscript for thesis defense, internal circulation, or a methods-oriented preprint. For an IEEE Transactions submission, complete the calibrated ANDES IBR benchmark and add full eigensolve comparisons for estimator accuracy, conversion ranking regret, and line-sensitivity finite differences.",
            "",
        ]
    )
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
