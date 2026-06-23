# Consolidated Paper Blueprint

Working title:

> A Damping-Region Spectral Estimator for Stability-Aware Planning of
> Inverter-Dominated Power Grids

Status: consolidated theory proposal. Claims are calibrated. Numerical claims
from random/small models are not final until IEEE 39-bus ANDES validation.

---

## 1. Core Research Question

How can the damping margin of an IBR-dominated grid with heterogeneous,
non-proportional damping be estimated cheaply and used to plan:

1. which synchronous generators can be converted to inverter-based resources,
2. where to add damping or virtual inertia cost-effectively,
3. which line reinforcements improve the actual damping margin rather than only
   increasing modal frequency?

The paper is about small-signal dynamic stability and planning. It is not an OPF
paper, not a transient-stability basin paper, and not a claim that graph
connectivity alone determines stability.

---

## 2. Model and Notation

Start from the linearized swing/network model:

```text
M theta_ddot + D theta_dot + L theta = p
```

where:

- `L` is the weighted active-power Laplacian / stiffness matrix,
- `M = diag(M_i)` is inertia or virtual inertia,
- `D = diag(D_i)` is damping or equivalent local damping.

Use inertial coordinates:

```text
L_tilde = M^(-1/2) L M^(-1/2)
D_tilde = M^(-1/2) D M^(-1/2)
```

Then:

```text
vartheta_ddot + D_tilde vartheta_dot + L_tilde vartheta = M^(-1/2) p
```

Let:

```text
L_tilde q_k = nu_k q_k
Q^T Q = I
Gamma = Q^T D_tilde Q
delta_k = Gamma_kk = q_k^T D_tilde q_k
```

The exact modal dynamics are governed by the QEP:

```text
P(s) v = (s^2 I + s D_tilde + L_tilde) v = 0
```

The exact damping margin is:

```text
zeta_min = min_j -Re(s_j) / |s_j|
```

over oscillatory complex roots.

---

## 3. Formal Theory

### Proposition 1: Proportional Damping Region

If `L_tilde` and `D_tilde` commute, then `Gamma` is diagonal and each mode
decouples:

```text
s^2 + delta_k s + nu_k = 0
```

For an underdamped mode:

```text
zeta_k = delta_k / (2 sqrt(nu_k))
```

Therefore, for a target margin `zeta_star`, the system satisfies the modal
margin requirement if every oscillatory mode lies in:

```text
R_diag(zeta_star) = {(nu, delta): delta >= 2 zeta_star sqrt(nu)}
```

This is the damping-region / damping-consensus picture.

Status:

- Verified mathematically.
- Known in spirit for homogeneous swing dynamics.
- Candidate contribution is the IBR planning geometry and its extension to
  heterogeneous modal damping, not the proportional formula itself.

### Proposition 2: First-Order Non-Proportional Approximation

For weak damping and well-separated modal frequencies, the diagonal estimate:

```text
zeta_diag = min_k delta_k / (2 sqrt(nu_k))
```

is the first-order approximation of the exact QEP margin.

The approximation can fail when:

```text
|Gamma_lk| is large relative to |sqrt(nu_l) - sqrt(nu_k)|
```

or when several low-damping modes form a frequency cluster.

This gives the diagnostic logic:

```text
if modal coupling is weak and modes are separated:
    use zeta_diag
else:
    solve reduced QEP on the coupled low-damping cluster
```

Status:

- Verified numerically in synthetic tests: diagonal is excellent for weak
  damping; it fails in near-degenerate/strongly coupled cases.
- Needs formal error bound or empirical diagnostic threshold before submission.

### Proposition 3: Reduced-QEP Correction

Let `I_r` be a selected set of modes containing:

- the lowest `zeta_diag` modes,
- any neighboring modes strongly coupled through `Gamma_lk`,
- any near-degenerate modes by `nu` or `sqrt(nu)`.

Let `V = Q[:, I_r]`. Define:

```text
L_r = V^T L_tilde V
D_r = V^T D_tilde V
```

Solve:

```text
(s^2 I_r + s D_r + L_r) w = 0
```

and estimate:

```text
zeta_adapt = min(zeta_diag, zeta_QEP_r)
```

or use `zeta_QEP_r` when the coupling diagnostic triggers.

Status:

- The reduced QEP is the correction, not always the base method.
- This is the clean version of the older C5 claim.

### Proposition 4: Damping-Aware Line Sensitivity

The old frequency sensitivity is:

```text
d nu_k / d w_e = (q_k,i / sqrt(M_i) - q_k,j / sqrt(M_j))^2
```

for edge `e=(i,j)`.

This only measures frequency/stiffness. It does not measure damping margin.
For the diagonal region:

```text
zeta_k = delta_k / (2 sqrt(nu_k))
```

so:

```text
d zeta_k / d w_e =
    (1 / (2 sqrt(nu_k))) d delta_k / d w_e
    - (delta_k / (4 nu_k^(3/2))) d nu_k / d w_e
```

Even if `d nu_k / d w_e > 0`, the second term is negative. Therefore, a line
reinforcement can increase frequency while decreasing damping ratio unless
`delta_k` increases enough.

For a simple modal eigenvalue:

```text
d delta_k / d w_e =
    q_k^T (d D_tilde / d w_e) q_k
    + 2 sum_{l != k}
        (q_l^T (d L_tilde / d w_e) q_k)
        (q_l^T D_tilde q_k)
        / (nu_k - nu_l)
```

If line reinforcement changes only `L_tilde`, then `d D_tilde / d w_e = 0`,
but `delta_k` can still change through eigenvector rotation. This is exactly
where non-proportional damping matters.

For the exact QEP, let `s` be the critical complex root and `v,y` the right and
left QEP eigenvectors:

```text
P(s) v = 0
y^* P(s) = 0
```

For parameter `rho`:

```text
ds / d rho =
  - y^* (s dD_tilde/drho + dL_tilde/drho) v
    / y^* (2s I + D_tilde) v
```

assuming `M` fixed and no `s^2` coefficient perturbation in normalized
coordinates. If `M` also changes, include `dL_tilde`, `dD_tilde`, and the
coefficient perturbation induced by the coordinate change.

For `s = sigma + j omega` and `ds/drho = a + j b`:

```text
d zeta / d rho = (sigma omega b - omega^2 a) / |s|^3
```

This is the analytic route for the damping-aware line map:

```text
S_zeta[e] = d zeta_min / d w_e
```

Status:

- Finite-difference verification is reported in the cross-audit.
- Analytic QEP sensitivity should be implemented and compared against finite
  differences before making it a formal contribution.

---

## 4. Paper Contributions

### C1. Damping-Region Framework

Claim:

A power-grid damping margin can be visualized and screened in the modal plane
`(nu_k, delta_k)` using the region:

```text
delta_k >= 2 zeta_star sqrt(nu_k)
```

This generalizes the intuition of consensus regions to inertial power-grid
damping margins.

What is new:

- Not the proportional formula itself.
- The contribution is the planning geometry for heterogeneous IBR damping and
  the diagnostic that tells when the diagonal region is insufficient.

Risk:

- Similar stability-region ideas are active in converter-dominated grid
  literature. Need targeted literature search before using "novel".

### C2. Adaptive Spectral Margin Estimator

Claim:

The damping margin can be estimated cheaply by:

1. computing diagonal modal damping over all modes,
2. diagnosing non-proportional modal coupling,
3. solving a reduced QEP only when coupling/near-degeneracy demands it.

Why it matters:

- It is cheaper and clearer than full eigensolve.
- It avoids overclaiming that reduced QEP is always necessary.
- It also avoids the false simplification that diagonal formulas always work.

Required validation:

- diagonal vs reduced-QEP vs full eigensolve,
- weak vs strong damping,
- separated vs clustered modes,
- synthetic grids and IEEE 39 ANDES.

### C3. Stability-Aware SG-to-IBR Conversion and Placement

Claim:

Planning actions should be ranked by post-action damping margin:

```text
score(action) = zeta_min_after_action / cost(action)
```

Actions include:

- converting one SG to GFL/GFM,
- adding virtual damping,
- adding virtual inertia,
- compensating a bus/controller parameter.

Why it matters:

- It answers a practical planning question.
- It avoids misleading proxies such as Fiedler participation or `lambda_2`.

Required validation:

- single conversion,
- multi-node progressive conversion,
- comparison against naive modal-participation rules,
- ANDES IEEE 39 validation.

### C4. Damping-Aware Mode-Line Reinforcement Map

Claim:

Line reinforcement should be ranked by:

```text
d zeta_min / d w_e
```

or by finite post-action `Delta zeta_min`, not by `d nu_k / d w_e` alone.

Mechanism:

Reinforcing a line generally increases modal frequency. Since:

```text
zeta = -sigma / sqrt(sigma^2 + omega^2)
```

increasing `omega` without increasing decay `-sigma` enough can lower `zeta`.
This is a damping-frequency tradeoff, not a topological Braess paradox.

Required validation:

- compare `d nu/dw` ranking vs `d zeta/dw` ranking,
- finite-difference check,
- analytic QEP-sensitivity check,
- full ANDES eigen-sensitivity check.

---

## 5. What Not to Claim

Do not claim:

- `lambda_2` or Fiedler mode is the damping margin,
- the highest-frequency mode is always critical under heterogeneous damping,
- line reinforcement is always stabilizing,
- the diagonal formula is exact under non-proportional damping,
- reduced QEP is always necessary,
- dynamic Braess paradox unless the mechanism is truly topological and not
  reduced modal damping,
- novelty before literature search.

---

## 6. Proposed Paper Structure

Target length: 10-14 pages conference/journal short paper; 20 pages for thesis
chapter.

### Abstract

Problem: IBR penetration makes damping heterogeneous and non-proportional.
Classical homogeneous formulas and graph-connectivity proxies are insufficient.

Contributions:

1. damping-region view,
2. adaptive diagonal-plus-QEP margin estimator,
3. SG-to-IBR / damping-placement planning rule,
4. damping-aware line reinforcement map.

Validation: synthetic tests plus IEEE 39-bus ANDES target.

### I. Introduction

- IBR transition and damping-margin problem.
- Why low inertia alone is not the full story: damping/control heterogeneity.
- Why `lambda_2`, Fiedler participation, and frequency-only reinforcement are
  insufficient.
- Contributions and honesty note: formulas are built on known modal/QEP theory;
  contribution is the calibrated estimator and planning use.

### II. Related Work

Organize by:

- low-inertia spectral performance metrics,
- consensus regions and MAS spectral conditions,
- non-proportional damping and QEP sensitivity,
- IBR small-signal stability/security regions,
- power-system damping/eigenvalue sensitivity.

Positioning:

- proportional damping formulas are known,
- stability/security regions are active,
- line/eigenvalue sensitivities are known,
- the candidate gap is the IBR damping-region estimator with adaptive QEP
  correction and damping-aware planning maps.

### III. Model and Damping-Region Theory

Contents:

- swing / small-signal model,
- inertial normalization,
- QEP,
- exact damping margin,
- Proposition 1: proportional damping region,
- Proposition 2: diagonal first-order approximation,
- coupling diagnostic.

Main figure:

- scatter of modes in `(sqrt(nu), delta)` or `(nu, delta)`,
- safe boundary `delta = 2 zeta_star sqrt(nu)`,
- critical modes highlighted.

### IV. Adaptive Spectral Estimator

Algorithm:

```text
Input: L, M, D, target zeta_star, r/base threshold
1. Build L_tilde, D_tilde.
2. Compute nu, Q, Gamma.
3. Compute zeta_diag for all modes.
4. Detect coupling / clustered low-damping modes.
5. If safe: return zeta_diag.
6. Else: solve reduced QEP on selected cluster.
7. Return margin and diagnostic label.
```

Experiments:

- proportional control,
- weak damping heterogeneity,
- strong heterogeneity,
- near-degenerate cases.

Report:

- median and p95 error,
- failure cases,
- when QEP correction is triggered.

### V. Planning Application: SG-to-IBR and Damping/Inertia Placement

Actions:

- convert candidate SG to GFL/GFM,
- add damping at node `i`,
- add virtual inertia at node `i`.

Scores:

```text
score_i = zeta_min(action_i)
cost_score_i = Delta zeta_min(action_i) / cost_i
```

Baselines:

- Fiedler participation,
- highest-frequency participation,
- inertia-only heuristics,
- homogeneous formula.

Validation:

- single action,
- greedy multi-action,
- ANDES 39-bus.

### VI. Damping-Aware Line Reinforcement

Start with the classical frequency map:

```text
d nu_k / d w_e
```

Then show why it is insufficient. Derive:

```text
d zeta_k / d w_e
```

and exact QEP sensitivity:

```text
ds / d w_e
```

Experiments:

- compare top lines by `d nu/dw` and `d zeta/dw`,
- show cases where `nu` rises but `zeta_min` falls,
- finite-difference validation,
- analytic sensitivity validation.

Main claim:

The corrected C4 is not "which line raises frequency"; it is "which line raises
the damping margin."

### VII. IEEE 39-Bus ANDES Validation

Required setup:

- IEEE 39 New England system,
- SG models: e.g., GENROU where available,
- IBR models: documented GFL/GFM/PLL models available in ANDES,
- penetration sweep,
- no hand-tuned gains unless explicitly labeled as sensitivity analysis.

Experiments:

1. margin estimator: full eigensolve vs diagonal vs adaptive QEP,
2. conversion ranking,
3. damping/inertia placement,
4. line reinforcement ranking,
5. ablations: frequency-only, Fiedler-only, homogeneous formula.

### VIII. Discussion

Include:

- what is exact,
- what is approximate,
- when the diagonal region fails,
- why QEP correction is needed,
- scalar Laplacian limitation under full dq/network dynamics,
- why G2/H-infinity and G3/switched topology are not central in this paper.

### IX. Conclusion

Restate:

- damping-region view,
- adaptive estimator,
- planning conversion/placement,
- damping-aware line map.

Future work:

- full dq/matrix-weighted Laplacian,
- H-infinity disturbance-shape margin,
- switched topology/dwell-time theory,
- delay-dependent control margins only if the project scope changes.

---

## 7. Figure and Table Plan

Figure 1:

- Conceptual damping region in `(nu, delta)` with safe/unsafe modes.

Figure 2:

- Diagonal vs adaptive QEP vs full eigensolve error across regimes.

Figure 3:

- Coupling diagnostic: off-diagonal modal damping / spectral gap predicts when
  diagonal fails.

Figure 4:

- SG-to-IBR conversion ranking: predicted vs true post-conversion margin.

Figure 5:

- Line reinforcement: frequency-only map vs damping-aware map.

Figure 6:

- IEEE 39 ANDES validation dashboard.

Table 1:

- Contribution status: exact, approximate, validated synthetic, validated ANDES.

Table 2:

- Baseline comparison: homogeneous formula, Fiedler, frequency-only, diagonal,
  adaptive QEP.

Table 3:

- Failure cases and what mechanism explains them.

---

## 8. Current Claim Labels

| Claim | Status |
| --- | --- |
| Proportional damping region | Verified, known in spirit |
| Diagonal all-mode estimator | Verified as first-order / weak-damping approximation |
| Adaptive QEP correction | Verified in synthetic tests, needs ANDES |
| SG-to-IBR ranking | Promising synthetic application, needs ANDES |
| Damping-aware line map | Mechanistically sound, finite-difference reported, analytic implementation pending |
| H-infinity hidden margin | Secondary, not central |
| Switched topology dwell-time | Future work / separate thesis |
| Matrix-weighted dq Laplacian | High-value future theory |

---

## 9. Immediate Implementation Checklist

1. Implement `estimate_margin_adaptive(L, M, D)`.
2. Implement modal coupling diagnostic.
3. Implement finite-difference `line_zeta_sensitivity`.
4. Implement analytic QEP eigenvalue sensitivity.
5. Validate analytic sensitivity against finite differences.
6. Re-run synthetic suites with fixed seeds and logged regimes.
7. Update master paper plan and abstract.
8. Build ANDES extraction pipeline.
9. Run IEEE 39 validation.
10. Run targeted literature search before any "novel" claim.

