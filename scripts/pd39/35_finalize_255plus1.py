from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def f(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def fmt(x, digits=9):
    return "NA" if not math.isfinite(x) else f"{x:.{digits}g}"


def finite(values):
    return [x for x in values if math.isfinite(x)]


def max_abs(rows, key):
    return max((abs(f(r, key)) for r in rows if math.isfinite(f(r, key))), default=math.nan)


def current_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unavailable"


audit = read_csv(RESULTS / "PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv")
ind = [r for r in audit if r["audit_layer"] == "independent" and r["status"] == "ok"]
tol = [r for r in audit if r["audit_layer"] == "tolerance" and r["status"] == "ok"]
tol_by_key = {(r["portfolio"], r["scenario"], r["tolerance"]): f(r, "alpha") for r in tol}
tol_diffs = []
for p, s, _ in {(r["portfolio"], r["scenario"], r["tolerance"]) for r in tol}:
    a10 = tol_by_key.get((p, s, "1.0e-10"), tol_by_key.get((p, s, "1e-10"), math.nan))
    a12 = tol_by_key.get((p, s, "1.0e-12"), tol_by_key.get((p, s, "1e-12"), math.nan))
    if math.isfinite(a10) and math.isfinite(a12):
        tol_diffs.append(abs(a10 - a12))
bad_ind = [r for r in ind if abs(f(r, "fd_alpha_difference")) > 1e-4]
sign_changes = sum((f(r, "reference_alpha") < 0) != (f(r, "independent_fd_alpha") < 0) for r in ind)
robust_changes = sum((f(r, "reference_alpha") <= -0.05) != (f(r, "independent_fd_alpha") <= -0.05) for r in ind)

tds = read_csv(RESULTS / "PD39_255PLUS1_HIGH_PLL_TDS.csv")
tds_ok = [r for r in tds if r["status"] == "ok"]
tds_groups = {}
for r in tds_ok:
    tds_groups.setdefault(r["case_id"], []).append(r)

holdout = read_csv(RESULTS / "PD39_255PLUS1_HOLDOUT_CENSUS.csv")
holdout_summary = read_csv(RESULTS / "PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv")
h0_counts = [int(r["h0_count"]) for r in sorted(holdout_summary, key=lambda r: r["condition"])]
h005_counts = [int(r["h005_count"]) for r in sorted(holdout_summary, key=lambda r: r["condition"])]
v8 = "30;32;33;34;35;36;37;38"
v8_h0 = [r["condition"] for r in holdout_summary if v8 in r["h0_portfolios"].split("|")]
v8_h005 = [r["condition"] for r in holdout_summary if v8 in r["h005_portfolios"].split("|")]

composition = read_csv(RESULTS / "PD39_255PLUS1_PENETRATION_COMPOSITION.csv")
composition_pairs = {}
for r in composition:
    composition_pairs.setdefault(r["pair_id"], []).append(r)
composition_effects = [f(r, "delta_alpha_a_minus_b") for r in composition if math.isfinite(f(r, "delta_alpha_a_minus_b"))]
composition_pair_summary = []
for pair_id, rows in composition_pairs.items():
    effects = finite([f(r, "delta_alpha_a_minus_b") for r in rows])
    composition_pair_summary.append({
        "pair_id": pair_id,
        "cardinality": int(rows[0]["cardinality"]),
        "alpha_source": rows[0]["alpha_source"],
        "min_delta_alpha": min(effects) if effects else None,
        "max_delta_alpha": max(effects) if effects else None,
        "max_abs_delta_alpha": max(map(abs, effects)) if effects else None,
    })

modal = read_csv(RESULTS / "PD39_255PLUS1_MODE_TRACKING.csv")
macs = finite([f(r, "mac_to_v8") for r in modal])
mode_families = sorted(set(r["critical_mode_family"] for r in modal if r["critical_mode_family"]))

baseline = read_csv(RESULTS / "PD39_CLASSICAL_BASELINES.csv")
seven = sorted((r for r in baseline if int(r["cardinality"]) == 7), key=lambda r: f(r, "m_9"))

closure = read_csv(RESULTS / "PD39_255PLUS1_CLOSURE_AUDIT.csv")
homotopy = read_csv(RESULTS / "PD39_255PLUS1_HOMOTOPY_AUDIT.csv")

gate_a = {
    "status": "FAIL",
    "audit_rows": len(audit),
    "tolerance_rows": len(tol),
    "independent_rows": len(ind),
    "tolerance_max_alpha_difference_1e10_vs_1e12": max(tol_diffs, default=math.nan),
    "independent_alpha_failures_over_1e-4": len(bad_ind),
    "independent_max_fd_alpha_difference": max_abs(ind, "fd_alpha_difference"),
    "descriptor_max_difference": max_abs(ind, "descriptor_alpha_difference"),
    "fd_descriptor_max_difference": max_abs(ind, "fd_descriptor_difference"),
    "sign_changes": sign_changes,
    "robust_classification_changes": robust_changes,
    "interpretation": "Tolerance and descriptor checks pass, but the frozen central-FD alpha check fails because the rightmost eigenvalue is highly conditioned/non-normal.",
}

gate_b = {
    "status": "PASS",
    "rows": len(tds),
    "ok_rows": len(tds_ok),
    "sign_consistent_rows": sum(r["sign_consistent"].lower() == "true" for r in tds_ok),
    "load_bus": int(tds_ok[0]["load_bus"]) if tds_ok else None,
    "pulse_start_s": 1.0,
    "pulse_end_s": 1.1,
    "groups": {
        key: {
            "mean_linear_alpha": mean(f(r, "linear_alpha") for r in rows),
            "mean_estimated_rate": mean(f(r, "estimated_rate_s_inv") for r in rows),
            "min_estimated_rate": min(f(r, "estimated_rate_s_inv") for r in rows),
            "max_estimated_rate": max(f(r, "estimated_rate_s_inv") for r in rows),
        }
        for key, rows in sorted(tds_groups.items())
    },
    "interpretation": "All 12 integrations completed; TDS preserves the V8 positive-growth versus 7/8 negative-decay ordering.",
}

gate_c = {
    "status": "FAIL_GENERALIZATION",
    "rows": len(holdout),
    "unique_portfolio_condition_keys": len({(r["portfolio"], r["condition"]) for r in holdout}),
    "portfolios": len({r["portfolio"] for r in holdout}),
    "conditions": len({r["condition"] for r in holdout}),
    "failed_equilibrium_or_spectrum_rows": sum(r["status"] != "ok" for r in holdout),
    "h0_counts_by_condition": h0_counts,
    "h005_counts_by_condition": h005_counts,
    "h0_total": sum(h0_counts),
    "h005_total": sum(h005_counts),
    "conditions_with_nonzero_h0": sum(x > 0 for x in h0_counts),
    "conditions_with_nonzero_h005": sum(x > 0 for x in h005_counts),
    "v8_h0_conditions": v8_h0,
    "v8_h005_conditions": v8_h005,
    "interpretation": "The exact discovery 255+1 structure does not generalize as a universal holdout identity; blockers appear in 12 of 24 conditions and V8 is not an H0.05 blocker in the holdout census.",
}

gate_d = {
    "status": "PARTIAL",
    "composition_rows": len(composition),
    "composition_finite_effects": len(composition_effects),
    "composition_max_abs_delta_alpha": max(map(abs, composition_effects)) if composition_effects else None,
    "composition_pair_summary": composition_pair_summary,
    "modal_rows": len(modal),
    "mode_families": mode_families,
    "mac_min": min(macs) if macs else None,
    "mac_max": max(macs) if macs else None,
    "homotopy_status": homotopy[0]["status"] if homotopy else "unavailable",
    "interpretation": "Aggregate-matched composition effects are present, while common-observable mode tracking is modal rather than a continuous homotopy; physical SG-to-GFL homotopy is blocked by discrete model semantics.",
}

gate_e = {
    "status": "BLOCKED",
    "items": len(closure),
    "all_statuses": sorted(set(r["status"] for r in closure)),
    "interpretation": "Exact K,D,Q objects, dimensions, units, determinant identity, and Q->-1 boundary are not defined by the installed model/API. No proxy was used.",
}

headline = {
    "campaign": "PD39 255+1 - final mechanism and generalization validation",
    "base_commit": "902403cf",
    "branch": "research/pd39-255plus1-mechanism-validation",
    "campaign_head_at_artifact_generation": current_commit(),
    "prereg_commit": "5f28df82",
    "final_case": "C",
    "gate_status": {"A": "FAIL", "B": "PASS", "C": "FAIL_GENERALIZATION", "D": "PARTIAL", "E": "BLOCKED"},
    "numerical_truth": gate_a,
    "high_pll_tds": gate_b,
    "holdout_census": gate_c,
    "penetration_composition_mode": gate_d,
    "network_closure": gate_e,
    "excluded_campaigns": ["weak nodes", "weak links", "structured-radius 128-direction", "repair optimization", "planners", "co-design", "IEEE-68", "EMT", "new controller search"],
    "ias_poster_redesign": "NO",
}

(RESULTS / "PD39_255PLUS1_HEADLINE.json").write_text(json.dumps(headline, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

claims = [
    ["H1", "A", "V8 is the only minimal 0.05 robustness blocker in frozen discovery", "PASS_DISCOVERY", "Frozen discovery has the exact 255+1 result; holdout generalization is separately false.", "results/pd39/portfolio_campaign/portfolio_scenario_results.csv; results/PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv"],
    ["H2", "A", "V8 true-instability cases are numerical-truth reproducible", "PARTIAL_FAIL", "Tolerance and descriptor checks agree, but 72/81 independent FD alpha comparisons exceed 1e-4; signs remain unchanged.", "results/PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv"],
    ["H3", "B", "High-PLL nonlinear TDS preserves fragile/robust ordering", "PASS", "12/12 cases complete and 12/12 sign-consistent; V8 mean rate is positive while both 7/8 controls decay.", "results/PD39_255PLUS1_HIGH_PLL_TDS.csv"],
    ["H4", "C", "255+1 structure generalizes across all 24 holdouts", "FAIL", "6144/6144 keys evaluated, but H0 total is 42 and H0.05 total is 47 across 12/24 conditions.", "results/PD39_255PLUS1_HOLDOUT_CENSUS.csv; results/PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv"],
    ["H5", "D", "Aggregate-matched portfolios retain composition-sensitive alpha differences", "PASS", "90 matched comparisons have finite alpha effects; aggregate MW/MVA differences are zero by construction.", "results/PD39_255PLUS1_PENETRATION_COMPOSITION.csv"],
    ["H6", "D", "7/8 to 8/8 is a physically trackable mode mechanism", "PARTIAL", "All 81 audited modes are classified electromechanical/control; common-observable MACs are portfolio-dependent, but continuous homotopy is blocked.", "results/PD39_255PLUS1_MODE_TRACKING.csv; results/PD39_255PLUS1_HOMOTOPY_AUDIT.csv"],
    ["H7", "E", "Exact network-closure identity is testable in installed model", "BLOCKED", "K,D,Q semantics and determinant/Q boundary are absent; no proxy constructed.", "results/PD39_255PLUS1_CLOSURE_AUDIT.csv"],
    ["K1", "A", "Tight tolerance settings are distinguishable from loose settings", "PASS", "243/243 tolerance rows completed; max 1e-10 versus 1e-12 alpha difference is below 1e-9.", "results/PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv"],
    ["K2", "B", "TDS event is implemented at the preregistered bus and window", "PASS", "Load bus 39, pulse 1.0-1.1 s, 12/12 integrations completed.", "src/pd39/confirmatory.jl; results/PD39_255PLUS1_HIGH_PLL_TDS.csv"],
    ["K3", "C", "Holdout census has exactly 6144 attempts", "PASS", "6144 rows and 6144 unique portfolio-condition keys across 256 portfolios and 24 conditions.", "results/PD39_255PLUS1_HOLDOUT_CENSUS.csv"],
    ["K4", "D", "Remaining MW/MVA/inertia are not treated as mechanism proof", "PASS", "Matched composition table reports aggregate metadata separately from alpha and mode fields.", "results/PD39_255PLUS1_PENETRATION_COMPOSITION.csv"],
    ["K5", "E", "No proxy closure construction is acceptable", "PASS", "Exact closure gate is explicitly BLOCKED and no Laplacian or arbitrary proxy is substituted.", "results/PD39_255PLUS1_CLOSURE_AUDIT.csv"],
]
with (RESULTS / "PD39_255PLUS1_MASTER_CLAIMS.csv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out)
    writer.writerow(["claim_id", "gate", "claim", "decision", "quantitative_evidence", "interpretation", "evidence_files"])
    writer.writerows(claims)

seven_lines = []
for r in seven:
    missing = next((b for b in [30, 32, 33, 34, 35, 36, 37, 38] if str(b) not in r["portfolio"].split(";")), "?")
    seven_lines.append(f"- missing {missing}: portfolio {r['portfolio']}; converted MW {r['converted_mw']}; m9 {fmt(f(r, 'm_9'), 12)}; nominal alpha {fmt(f(r, 'nominal_alpha'), 12)}")

tds_lines = []
for key, rows in sorted(tds_groups.items()):
    tds_lines.append(f"- {key}: mean linear alpha {fmt(mean(f(r, 'linear_alpha') for r in rows), 10)}; mean nonlinear estimated rate {fmt(mean(f(r, 'estimated_rate_s_inv') for r in rows), 10)}")

holdout_lines = [f"{r['condition']}: H0={r['h0_count']}, H0.05={r['h005_count']}, failed={r['failed']}" for r in sorted(holdout_summary, key=lambda r: r['condition'])]

report = f"""# PD39 255+1 - Final mechanism and generalization validation

Campaign base: commit 902403cf
Branch: research/pd39-255plus1-mechanism-validation
Preregistration commit: 5f28df82

FINAL CASE: C

## 1. Executive verdict

This campaign did execute the requested new four-gate validation. It does not validate the strongest version of the 255+1 story. The frozen discovery pattern remains exact in the original nine conditions, and the nonlinear high-PLL TDS preserves the fragile/robust ordering. However, the preregistered independent central-FD alpha check fails for {len(bad_ind)}/81 cases; the complete 24-condition census produces {sum(h0_counts)} H0 blockers and {sum(h005_counts)} H0.05 blockers; the physical SG-to-GFL homotopy is unavailable; and the exact K,D,Q closure test is BLOCKED by missing model semantics. Under the preregistered decision rule this is FINAL CASE C: the claims must be narrowed and the IAS poster should not be redesigned around a validated universal mechanism.

## 2. Gate A - numerical truth audit

The audit contains 9 portfolios x 9 scenarios x 3 tolerances = 243 tolerance rows and 81 independent rows. All 243 tolerance rows completed. The maximum alpha difference between 1e-10 and 1e-12 was {fmt(gate_a['tolerance_max_alpha_difference_1e10_vs_1e12'], 12)} s^-1, with no sign or 0.05-classification changes in that tolerance comparison. Descriptor generalized eigenvalues also agree with the reference reduction; the maximum descriptor difference was {fmt(gate_a['descriptor_max_difference'], 12)} s^-1.

The preregistered independent central finite-difference Jacobian was not a successful independent alpha validation: {len(bad_ind)}/81 comparisons exceeded 1e-4 s^-1, with maximum difference {fmt(gate_a['independent_max_fd_alpha_difference'], 12)} s^-1. There were {sign_changes} true-stability sign changes but {robust_changes} engineering-margin classification changes. A one-case step probe showed that the full Jacobian norm can be close while the rightmost eigenvalue is highly sensitive; the V8 nominal eigenvector condition was approximately 5.8e4. This explains the numerical behavior but does not erase the preregistered failure.

Gate A: FAIL on the frozen independent-alpha criterion; tolerance and descriptor subchecks PASS.

## 3. Gate B - high-PLL nonlinear TDS

The corrected implementation applies the common +1% P/Q constant-power-factor pulse at deterministic load bus 39 from 1.0 to 1.1 s and integrates to 20 s. Exactly 12 simulations completed, all with sign consistency.

{chr(10).join(tds_lines)}

V8 has positive estimated growth in all four high-PLL cases. Both 7/8 controls decay, with the best 7/8 more strongly damped than the worst 7/8. Gate B: PASS for sign and relative-order consistency. The earlier 50 ms implementation is not used as a result.

## 4. Gate C - complete holdout census

The census contains exactly {len(holdout)} rows, {len({(r['portfolio'], r['condition']) for r in holdout})} unique portfolio-condition keys, 256 portfolios, and 24 conditions. There are {sum(r['status'] != 'ok' for r in holdout)} failed equilibrium/spectrum attempts; all are retained and are not silently scored as instability.

Condition-wise counts:

{chr(10).join('- ' + x for x in holdout_lines)}

Totals are H0={sum(h0_counts)} and H0.05={sum(h005_counts)}. Nonzero blockers occur in 12/24 conditions for both definitions. V8 appears in H0 only at {', '.join(v8_h0) if v8_h0 else 'none'} and appears in H0.05 at {', '.join(v8_h005) if v8_h005 else 'none'}; therefore the exact V8 H0.05 identity does not generalize. Gate C: FAIL_GENERALIZATION, despite the census execution itself passing its row-count integrity checks.

## 5. Gate D - penetration, composition, and mode mechanism

The eight 7/8 predecessors have the following frozen discovery robust margins:

{chr(10).join(seven_lines)}

The 6/8 and 7/8 matched comparisons are aggregate-matched: converted MW and IBR MVA differences are zero within each pair, while alpha differences remain nonzero and scenario-dependent. The 90-row table distinguishes frozen-discovery alpha for 6/8 pairs from AD-reference alpha for the newly audited 7/8 pairs; no finite-difference alpha is used for this interpretation.

All 81 newly tracked critical modes are classified as electromechanical/control by the installed modal participation rule. V8 is a low-frequency oscillatory critical mode at roughly 0.32-0.37 Hz in the four high-PLL cases, with positive alpha roughly 0.110-0.114 s^-1. The 7/8 critical modes are generally distinct in the common bus-voltage observable: MAC values range from {fmt(min(macs), 6)} to {fmt(max(macs), 6)} where available, and several 7/8 critical modes are real/near-zero-frequency or have different oscillatory frequencies. The evidence supports a physical modal distinction, but does not establish a continuous SG-to-GFL homotopy or a single continuously tracked pole from 7/8 to 8/8.

Gate D: PARTIAL. Composition effects and modal tracking are supported; the physical homotopy portion is BLOCKED by discrete model semantics.

## 6. Gate E - exact network closure

The installed repository/API does not expose exact K, D, and Q objects with documented dimensions, units, construction path, determinant identity, or Q -> -1 boundary semantics. The audit therefore marks all five requested closure items BLOCKED. No graph Laplacian, arbitrary proxy, or fabricated identity was substituted.

Gate E: BLOCKED.

## 7. Claims and interpretation

The defensible result is narrower than a universal mechanism claim:

1. The original frozen discovery set contains an exact 255+1 robustness blocker pattern under its stated conditions.
2. Four high-PLL nonlinear TDS cases support the ordering V8 fragile, worst 7/8 intermediate, best 7/8 robust.
3. Aggregate-matched comparisons show composition-sensitive alpha differences, so converted MW/MVA alone do not explain every comparison.
4. The holdout census rejects the statement that the same H0 or H0.05 identity persists across all 24 designed conditions.
5. Independent finite-difference alpha validation is numerically unresolved/failing under the frozen step protocol because of strong eigenvalue sensitivity; this prevents a clean numerical-truth headline.
6. Exact closure theory and a physical continuous SG-to-GFL homotopy remain unavailable in this installed model.

The work that was intentionally not run is weak-node analysis, weak-link analysis, the 128-direction structured-radius campaign, repair optimization, planners, co-design, IEEE-68, EMT, and new controller search.

## 8. Final decision

Numerical audit: FAIL on the preregistered independent-alpha gate.
High-PLL TDS: PASS.
H0 over 24 holdouts: {sum(h0_counts)} total; counts [{', '.join(map(str, h0_counts))}].
H0.05 over 24 holdouts: {sum(h005_counts)} total; counts [{', '.join(map(str, h005_counts))}].
Penetration/composition/mixed verdict: composition-sensitive effects are present in aggregate-matched comparisons; no mixed intervention was run or claimed.
Mode mechanism: physical electromechanical/control modal distinction is supported, but continuous 7/8 -> 8/8 homotopy is blocked and the independent FD alpha audit is unresolved.
Closure: BLOCKED.
IAS poster redesign: NO.

FINAL CASE: C
"""

(DOCS / "PD39_255PLUS1_FINAL_REPORT.md").write_text(report, encoding="utf-8")

handoff = f"""# PD39 255+1 ChatGPT handoff

La campaña nueva se ejecutó en la rama `research/pd39-255plus1-mechanism-validation`, creada desde `902403cf`. No se usó ni se modificó la rama confirmatoria antigua.

## Resultado ejecutivo

`FINAL CASE: C`

La lectura correcta es parcial y negativa en puntos esenciales: Gate B pasa, pero Gate A falla su criterio independiente de alpha; Gate C no generaliza el patrón exacto; Gate D es parcial; Gate E está bloqueado por semántica/API ausente. Por eso la recomendación para el poster IAS es **NO** rediseñarlo alrededor de una afirmación universal 255+1.

## Gates

- Numerical audit: FAIL. 243/243 filas de tolerancia y 81/81 auditorías independientes se ejecutaron. La comparación 1e-10 vs 1e-12 fue estable, y el descriptor generalizado concordó; pero {len(bad_ind)}/81 diferencias de alpha por Jacobiano FD excedieron 1e-4 s^-1 (máximo {fmt(gate_a['independent_max_fd_alpha_difference'], 12)}). No hubo cambios de signo, pero sí {robust_changes} cambios de clase 0.05.
- High-PLL TDS: PASS. 12/12 casos, bus 39, pulso +1% P/Q de 1.0 a 1.1 s, integración a 20 s, 12/12 con signo consistente. Medias de tasa: V8 {fmt(gate_b['groups']['V8']['mean_estimated_rate'], 10)}, best 7/8 {fmt(gate_b['groups']['best_7of8']['mean_estimated_rate'], 10)}, worst 7/8 {fmt(gate_b['groups']['worst_7of8']['mean_estimated_rate'], 10)} s^-1.
- H0 over 24 holdouts: total {sum(h0_counts)}; vector [{', '.join(map(str, h0_counts))}].
- H0.05 over 24 holdouts: total {sum(h005_counts)}; vector [{', '.join(map(str, h005_counts))}].
- Penetration/composition/mixed: 90 comparaciones aggregate-matched con diferencias MW/MVA agregadas nulas y efectos alpha no nulos; no se ejecutó intervención mixta y no se reclama una ventaja mixta.
- Mode mechanism: critical modes classified electromechanical/control. V8 high-PLL is oscillatory at about 0.32-0.37 Hz with positive alpha; 7/8 modes are portfolio-dependent and often distinct in common bus-voltage MAC. Homotopy continua SG->GFL: BLOCKED.
- Closure: BLOCKED. No existen objetos K,D,Q exactos documentados ni identidad determinant/Q->-1 utilizable; no se fabricó proxy.
- IAS poster redesign: NO.

## Los ocho 7/8 margins de discovery

{chr(10).join(seven_lines)}

## Artefactos

- `docs/PD39_255PLUS1_FINAL_REPORT.pdf`
- `docs/PD39_255PLUS1_FINAL_REPORT.md`
- `results/PD39_255PLUS1_HEADLINE.json`
- `results/PD39_255PLUS1_MASTER_CLAIMS.csv`
- `results/PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv`
- `results/PD39_255PLUS1_HIGH_PLL_TDS.csv`
- `results/PD39_255PLUS1_HOLDOUT_CENSUS.csv`
- `results/PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv`
- `results/PD39_255PLUS1_PENETRATION_COMPOSITION.csv`
- `results/PD39_255PLUS1_MODE_TRACKING.csv`
- `results/PD39_255PLUS1_HOMOTOPY_AUDIT.csv`
- `results/PD39_255PLUS1_CLOSURE_AUDIT.csv`
- `deliverables/PD39_255PLUS1_CHATGPT_UPLOAD.zip`

No push was performed. The excluded campaigns (weak nodes/links, structured radius, repairs, planners, co-design, IEEE-68, EMT, new controller search) were not run.

FINAL CASE: C
"""
(DOCS / "PD39_255PLUS1_CHATGPT_HANDOFF.md").write_text(handoff, encoding="utf-8")

print("wrote final report markdown, handoff, headline JSON, and master claims CSV")
