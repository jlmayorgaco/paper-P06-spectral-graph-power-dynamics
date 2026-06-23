# 24-Page Paper Refinement Strategy

Date: 2026-06-08

Status: planning and claim calibration after reviewing
`C:\Users\walla\Downloads\Paper_Final_RationalGraphFilter.pdf` and the current
IEEE draft in `paper_ieee_transactions/main.tex`.

## 1. Review of Claude's PDF

The PDF is not a complete paper. It is a compact four-page concept note. Its
value is that it gives the right unifying story:

1. IBR damping is non-proportional because inverter control contributes damping
   and memory not aligned with inertia.
2. The second-order pole correction separates a network term from a
   control-coupling term.
3. The same correction can be written as a rational graph filter of the pair
   `(L_tilde, D_tilde)`.
4. Weak nodes, weak links, and Braess-like topology reversals should be
   sensitivity-based applications, not the core novelty.
5. IEEE 39-bus ANDES validation is required before final benchmark claims.

The current `paper_ieee_transactions/main.tex` already contains most of the
honest estimator and validation logic. What it lacks for a long 24-page version
is:

- a formal graph-filter section with definitions and proof,
- a step-by-step derivation from DAE/power-system model to QEP surrogate,
- a clean separation between known theory and original formulation,
- an explicit comparison matrix against classical and recent methods,
- a reproducible IEEE 39 validation harness,
- a longer results plan with falsification criteria and tables ready to fill.

## 2. Best Paper Option

The strongest paper is not "weak nodes" and not "Braess." Those are
applications. The strongest paper is:

> Separating Network and Inverter-Control Effects on the Damping Margin of IBR
> Grids: A Rational Graph-Filter Perturbation Approach

The central claim should be:

> The damping margin of a non-proportionally damped IBR grid can be decomposed
> into a zero-order network prediction and a second-order inverter-control
> correction. That correction is a rational graph filter of the dynamic operator
> pair `(L_tilde, D_tilde)`, and it gives a diagnostic workflow for when a cheap
> spectral screen is trustworthy and when a QEP fallback is required.

This is a defensible power-systems contribution because the mathematical tools
are known, but the formulation, interpretation, and validation target are
specific to IBR damping-margin planning.

## 3. Contribution Hierarchy

### Central Contribution: second-order network-vs-control decomposition

Verified in reduced models:

```text
zeta_min ~= zeta_net^(0) + Delta zeta_ctrl^(2)
```

The pole correction is

```text
Delta s_c =
  s0^2 / (2 s0 + delta_c)
  sum_{ell != c} E_{c ell} E_{ell c}
  / (s0^2 + delta_ell s0 + nu_ell).
```

Claim status:

- verified in surrogate QEP ensembles,
- verified in the reduced seven-state PLL model,
- not yet verified in IEEE 39 ANDES,
- standard perturbation theory; originality is in IBR damping-margin use and
  network-vs-control interpretation.

### Second Contribution: rational graph-filter form

Define

```text
H(Lambda; s0) = diag(1 / (s0^2 + delta_ell s0 + nu_ell)).
```

Then

```text
Delta s_c =
  s0^2 / (2 s0 + delta_c) [E H(Lambda; s0) E]_{cc}.
```

This is a rational filter of the operator pair `(L_tilde, D_tilde)`.

Claim status:

- exact algebraic reformulation of the second-order correction,
- useful for graph signal processing framing,
- not a separate mechanism and should not be oversold.

### Third Contribution: adaptive estimator with trigger and QEP fallback

Workflow:

```text
diagonal all-mode screen
-> second-order correction if coupling/degeneracy diagnostic is safe
-> reduced QEP if the diagnostic rejects the perturbation
-> full eigensolve for benchmark reference
```

Claim status:

- verified as a useful structure in synthetic tests,
- trigger is conservative and empirical, not yet a tight theorem,
- QEP fallback is mandatory because second order can worsen tail cases.

### Fourth Contribution: weak nodes, weak links, and Braess-like diagnostics

Weak nodes:

```text
S_i = |partial zeta_min / partial rho_i|
```

where `rho_i` can be local damping, virtual inertia, droop, or PLL gain.

Weak links:

```text
S_ij = |partial zeta_min / partial w_ij|.
```

Braess-like damping reversal:

```text
Delta nu_c > 0   and   Delta zeta_min < 0.
```

Claim status:

- sensitivity definitions are correct,
- QEP sensitivity is known prior art,
- application to damping-aware line/topology decisions in IBR grids is
  secondary and must be compared against frequency-only and DRI/FDI baselines,
- do not call this a new Braess paradox unless a topological feedback mechanism
  is proved.

## 4. What Is Known Theory vs. Original Development

Known theory:

- proportional damping modal formula,
- inertial Laplacian modal decomposition,
- QEP formulation of second-order dynamics,
- eigenvalue and damping-ratio sensitivity,
- non-proportional damping modal perturbation,
- graph polynomial filters of a Laplacian,
- Braess-like effects in power-grid topology.

Original or plausibly original in this project:

- applying the closed-form second-order correction to the critical
  damping-margin pole of IBR grids,
- interpreting the zero-order term as the network margin and the second-order
  term as inverter-control coupling cost,
- expressing that correction as a rational graph filter of `(L_tilde,D_tilde)`
  and using it to explain why polynomial filters of `L_tilde` are insufficient,
- embedding the correction in a diagnostic estimator with QEP fallback,
- using the same margin object to define weak nodes, weak links, and
  damping-aware topology reversals.

## 5. Proposed 24-Page Structure

### Page 1: Title, abstract, contribution boundary

State the paper as a method and validation study. Avoid "new linear algebra."

### Pages 2-3: Introduction

Motivation:

- IBRs reduce inertia and introduce non-proportional damping.
- Frequency-strengthening actions can reduce damping margin.
- Planning needs a cheap but self-aware margin estimator.

Research question:

```text
How much damping margin comes from the network, how much is moved by inverter
control, and can this split support planning diagnostics?
```

### Pages 4-5: Related work and prior-art boundary

Organize by baseline category:

1. low-inertia/proportional damping spectral formulas,
2. GFM/GFL placement by damping or composite indices,
3. QEP and eigenvalue sensitivity methods,
4. non-proportional damping perturbation in structural dynamics,
5. graph signal processing and graph filters,
6. Braess/topology effects in power systems.

Explicitly state what is not claimed as new.

### Pages 6-8: Power-system model to second-order surrogate

Start from the linearized differential-algebraic model:

```text
E xdot = f_x x + f_y y
0      = g_x x + g_y y
```

Eliminate algebraic variables:

```text
xdot = A_sys x.
```

Then state the reduced second-order surrogate:

```text
M theta_ddot + D theta_dot + L theta = p.
```

Define:

```text
L = partial P_e / partial theta
M = inertia/virtual inertia
D = local damping/droop/equivalent PLL damping when physically meaningful
```

Be explicit that the surrogate is an approximation to be checked against full
ANDES poles.

### Pages 9-10: Damping region and diagonal all-mode screen

Derive:

```text
vartheta = M^(1/2) theta
L_tilde = M^(-1/2) L M^(-1/2)
D_tilde = M^(-1/2) D M^(-1/2)
```

Then:

```text
Gamma = Q^T D_tilde Q = Delta + E.
```

If `E=0`:

```text
zeta_k = delta_k / (2 sqrt(nu_k)).
```

Important: use all modes, not only Fiedler.

### Pages 11-13: Original theory, second-order correction

Derive the Schur-complement expansion step by step:

```text
P(s) = s^2 I + s(Delta + E) + Lambda.
```

For critical mode `c`, eliminate the other modal coordinates:

```text
p_c(s) - B(s) R(s)^(-1) C(s) = 0.
```

Expand about `s0`:

```text
p_c(s0) = 0
p_c'(s0) Delta s_c + second-order coupling = 0.
```

Obtain the correction formula and define corrected damping ratio.

### Pages 14-15: Rational graph-filter interpretation

Define polynomial graph filters and prove insufficiency:

```text
h(L_tilde) q_k = h(nu_k) q_k.
```

Therefore they cannot create modal coupling if the missing information is
`E = offdiag(Q^T D_tilde Q)`.

Then define the rational pair-filter:

```text
Delta s_c = alpha(s0) [E H(Lambda;s0) E]_{cc}.
```

Explain the physics:

- `E` injects non-proportional damping coupling,
- `H` weights coupling by modal dynamic separation,
- the critical pole selects the relevant frequency scale.

### Pages 16-17: Adaptive estimator and validity diagnostics

Algorithm:

1. compute `L_tilde` eigenpairs,
2. compute `Gamma`,
3. compute diagonal all-mode screen,
4. compute second-order correction,
5. compute coupling/degeneracy indicator,
6. choose diagonal, second-order, reduced QEP, or full eigensolve.

Report failure modes:

- clustered modes,
- strong off-diagonal damping,
- pole switching,
- surrogate mismatch.

### Pages 18-19: Weak nodes and weak links

Define intervention scores:

```text
S_i = |partial zeta_min / partial rho_i|
S_e = |partial zeta_min / partial w_e|.
```

Use QEP sensitivity:

```text
ds/d rho =
 - y^* (s dD/drho + dL/drho) v
 / y^* (2sI + D) v.
```

Then:

```text
d zeta/d rho =
 (sigma omega b - omega^2 a) / |s|^3
```

for `s=sigma+i omega`, `ds/d rho=a+i b`.

### Page 20: Braess-like damping topology reversal

State the diagnostic condition:

```text
Delta nu_c > 0 or Delta lambda_2 > 0,
Delta zeta_min < 0.
```

Frame honestly:

- not a new Braess paradox by itself,
- a damping-ratio reversal under topology/stiffness action,
- caused by decay-frequency geometry and non-proportional damping.

### Pages 21-23: IEEE 39 validation and comparisons

Required experiments:

1. estimator accuracy vs. full ANDES eigenanalysis,
2. network-vs-control decomposition vs. IBR penetration,
3. trigger precision/recall,
4. weak-node ranking regret,
5. weak-link ranking regret,
6. Braess-like reversal count and mechanism audit,
7. runtime/computational cost.

Baselines:

- homogeneous damping formula,
- Fiedler-only estimate,
- diagonal all-mode estimate,
- polynomial graph filters of `L_tilde`,
- low-frequency reduced QEP,
- damping-selected reduced QEP,
- full QEP/eigensolve,
- DRI/FDI GFM placement index,
- finite-difference full-eigensolve sensitivity,
- published QEP sensitivity formulas.

### Page 24: Discussion, limitations, and claim table

Include a table:

```text
Claim | Verified now | Not yet proven | Required IEEE 39 evidence
```

Include falsification criteria:

- if surrogate does not track ANDES poles, claim reduced-model screening only,
- if second order does not beat diagonal, claim interpretive decomposition,
- if weak-link rankings match frequency rankings, demote C4 to cautionary note,
- if trigger has false-safe cases, require full QEP fallback.

## 6. IEEE 39 Data Products Required

The final paper should produce:

1. `ieee39_case_metadata.json`: buses, generators, IBR replacements, controller
   gains, penetration levels.
2. `estimator_table.csv`: true pole, true margin, diagonal, second-order, QEP,
   trigger label, relative error.
3. `decomposition_table.csv`: network term, control correction, full margin,
   penetration level.
4. `weak_node_table.csv`: candidate bus, sensitivity, full post-action margin,
   regret.
5. `weak_link_table.csv`: line, frequency sensitivity, damping sensitivity,
   finite-difference margin change, reversal flag.
6. `braess_audit_table.csv`: action, `Delta nu`, `Delta lambda2`,
   `Delta zeta`, pole switching flag, accepted/rejected mechanism.
7. Figures:
   - damping-region scatter,
   - estimator error bars,
   - network-vs-control decomposition vs. penetration,
   - weak-node map,
   - weak-link map,
   - frequency-vs-damping reinforcement scatter,
   - trigger confusion matrix.

## 7. Decision Logic for Submission

For IEEE Transactions:

- require IEEE 39 ANDES validation,
- require comparisons against at least five baselines,
- require all claims labeled by evidence level,
- require reproducible scripts and case metadata.

For a master-thesis chapter or preprint:

- the current reduced-model evidence is enough to write the full theory,
- IEEE 39 can be a dedicated validation chapter,
- weak nodes/links/Braess can be included as proposed applications with partial
  evidence.

Recommended route:

1. Build the 24-page thesis/preprint version first.
2. Run IEEE 39 validation.
3. Split into one Transactions paper centered on the rational-filter estimator
   and one planning paper if weak nodes/links/Braess results are strong.

