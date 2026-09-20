# TX4 Robustness Audit Handoff

## Git
branch: research/tx4-robustness-statistics-audit  
parent: f64db0004026ceafdb08dd13b5e2ff59d6060742  
HEAD: d0764901065f11fe8e4bc565a5ec7167b863aec3  
commits: 04a3230a, cd42a7fb, d0764901  
runtime: exact fallback checkpoint wall time 6723.4 s including resume  
no push: YES

## Original problems found
- H4_PRESENT bug? YES. Legacy H4_PRESENT tested H4 alpha_EM instability rather than H4 membership in H0.
- NONCOMPOSABLE bug? YES. Legacy NONCOMPOSABLE represented solver failure instead of composite minimal blockers.
- eta_H4 issue? YES. Legacy eta was a proper-subset alpha margin, not the registered local/collective sigma-minimum diagnostic.
- surrogate use? YES. 302,460 proper-subset rows were surrogate-tier and failed the strict minimality gate.

## Data provenance
total rows: 329,440 legacy master rows  
exact rows: 145,536 exact fallback rows, plus 26,980 exact rows retained in the legacy master  
surrogate rows: 302,460  
all-16 exact conditions: 9,096 exact fallback conditions, plus 426 legacy calibration conditions

## Exact H4 instability
QMC count / N / coverage: 2,367 / 4,096 / 0.577881  
MC count / N / probability: 2,907 / 5,000 / 0.581400, Wilson 95% CI [0.567668, 0.595007]

## Surrogate validation
false minimal: 18  
missed minimal: 14  
antichain exact match: 0.701878  
EXACT_H4 agreement: 0.924883

## Was exact QMC+MC recomputation required?
YES. The strict surrogate gate failed.

Exact portfolio evaluations: 145,536

## Correct minimality results

QMC: H4_PRESENT=0.113770; EXACT_H4=0.113770; NONCOMPOSABLE=0.427246. KAPPA PMF is in `results/TX4_QMC_EXACT_BLOCKER_SUMMARY.csv`.

MC: H4_PRESENT=0.115400, Wilson CI [0.106838, 0.124553]; EXACT_H4=0.115400, Wilson CI [0.106838, 0.124553]; NONCOMPOSABLE=0.436600, Wilson CI [0.422907, 0.450391]. KAPPA PMF is in `results/TX4_MC_EXACT_BLOCKER_SUMMARY.csv`.

## delta_H4

QMC median=-0.085269 s^-1; MC median=-0.084810 s^-1 with BCa 95% CI [-0.088820, -0.081605].

## eta/local/collective

The exact stratified stratum contains 128 rows, zero failures, median eta=0.111537, minimum local sigma=0.294663, and median collective sigma=0.111537. All eta values are nonnegative.

## Sensitivity provenance

Morris remains approximate screening from the retained campaign. Sobol remains surrogate-based. Neither is used as exact primary minimality evidence.

## g_star provenance

g_star is retained as an exact-H4 local phase-boundary diagnostic, not a physical robust radius.

## Invariant checks

All pass: YES. The machine-readable file is `results/TX4_ROBUSTNESS_INVARIANT_CHECKS.csv` with 32 PASS rows.

## FINAL CASE

CASE B: H4 instability remains a strong exact benchmark result, while corrected exact minimality is substantially smaller and composite blockers are common.

## Strongest corrected result

The frozen reduced DAE supports exact H4 instability coverage of 57.8% QMC and 58.1% MC, but exact H4 minimality is 11.4% in both designs and NONCOMPOSABLE occurs in 42.7% QMC and 43.7% MC.

## Strongest previous claim that was wrong

The legacy H4_PRESENT percentages were presented as H4 presence/minimality even though they measured H4 instability; the legacy zero NONCOMPOSABLE result also used the wrong semantic definition.

## Exact wording allowed for poster

"In the frozen reduced IEEE-39 matched-policy DAE, exact QMC/MC screening found H4 instability in 57.8%/58.1% of the declared bounded design, but exact inclusion-minimal H4 blockers in 11.4%/11.5%; composite minimal blockers occurred in 42.7%/43.7%."

## Exact wording forbidden

Do not call these physical population probabilities, universal IEEE-39 claims, EMT or hardware certificates, or a physical robust radius. Do not reuse the legacy 57.8%/58.1% values as H4_PRESENT minimality.

## Paper implication

The corrected paper must separate H4 instability from minimal-blocker presence, state the exact fallback and provenance, report composite blockers and KAPPA, and retain the surrogate and bounded-box limitations.
