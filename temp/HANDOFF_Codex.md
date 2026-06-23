# Handoff: Spectral Damping-Margin Analysis for IBR Power Networks

This file is the first thing Codex should read before continuing this project.

Project goal: find a real, honest, technically defensible contribution at
master's thesis level / publishable-paper level on inverter-based-resource
(IBR) power-system stability using spectral graph theory, algebra, and
networked dynamical systems.

The priority is rigor, not impressive-looking results.

---

## 1. Honesty Protocol

Follow these rules without exception.

### Verify the mechanism, not the number

A large value, low error, high gain, or high percentage does not confirm the
hypothesis that motivated it. Before claiming that something works, show why it
works and rule out plausible alternative causes.

Most false positives in this project turned out to be low-damping effects or
modeling artifacts, not the hypothesized phenomenon.

### Compare against the correct control

When comparing methods or conditions, equalize everything else that matters:
damping ratio, system size, perturbation scale, operating point, model fidelity,
and metric definition.

Do not compare apples to oranges and then claim a mechanism.

### Separate verified from original

Numerically correct does not mean novel. Label every result as one of:

- verified and original,
- verified but already known / published,
- not verified.

Never use "original", "new", or "novel" unless a literature search has been
done for that exact claim. If the search has not been done, write:

> I have not confirmed that this is original.

### Do not revive dead hypotheses

Section 6 lists the hypotheses already refuted. If current reasoning points
toward one of them, stop and check why it failed before continuing.

### Report negative results

If a hypothesis does not hold, say so clearly and explain why. A verified
negative result is more valuable than a doubtful positive result.

Do not reframe a failed hypothesis as a partial success just to make it sound
better.

### Separate what is measured from what the model permits

Numbers from simplified models, 3-node examples, or toy PLL models are
illustrative only. They are not quantitative evidence for a paper.

Nothing is conclusive until validated on IEEE 39-bus with realistic ANDES
models for SG, GFL, GFM, PLL, and network dynamics as appropriate.

### Use calibrated epistemology

When uncertain, say so. Distinguish:

- what was tested,
- what is assumed,
- what is expected,
- what remains unknown.

For every reported result, include:

- what was verified,
- what it was compared against,
- what was not tested,
- confidence level.

If any of those four items are missing, the result is incomplete.

---

## 2. Current State

After the latest cross-audit, the paper should be framed around four calibrated
contributions. The older statement "the reduced QEP is always the principal
contribution" is too strong. The correct story is an adaptive spectral margin
framework:

- C1, theory frame: a damping-region / damping-consensus view. In the
  inertial-Laplacian modal plane `(nu_k, delta_k)`, the diagonal modal estimate
  satisfies `zeta_k ~= delta_k / (2 sqrt(nu_k))`, where
  `delta_k = q_k^T D_tilde q_k`. The safe region for a target margin
  `zeta_star` is `delta_k >= 2 zeta_star sqrt(nu_k)`. This is exact in the
  proportional/commuting case and a first-order diagnostic in weak, separated
  non-proportional cases. It is not a universal exact characterization of a
  non-proportional QEP.

- C2, main method: an adaptive diagonal-plus-reduced-QEP estimator. Use the
  diagonal estimate over all modes as the cheap baseline; diagnose modal
  coupling using the off-diagonal entries of `Q^T D_tilde Q` and spectral
  gaps; solve a damping-selected reduced QEP only for clustered/strongly coupled
  low-damping modes. This reconciles two facts: the diagonal formula can be very
  accurate for weak damping, while the reduced QEP is needed in degenerate or
  strongly non-proportional cases.

- C3, planning application: use the estimator to rank SG-to-IBR conversions and
  inertia/damping placement by post-action `zeta_min`, not by Fiedler
  participation or `lambda_2`. Single-node conversion looked promising in
  synthetic tests; multi-node progressive conversion and ANDES validation remain
  untested.

- C4, corrected line-reinforcement map: the old C7 frequency map
  `d nu_k / d w_ij` is useful but incomplete. Reinforcing a line can raise
  modal frequency while lowering damping ratio if modal damping does not rise
  enough. The paper should rank lines by `d zeta_min / d w_ij` or by finite
  post-action `Delta zeta_min`, with `d nu_k / d w_ij` reported only as a
  frequency component of the mechanism.

Main bottleneck: none of these claims is conclusive until validated on a real
benchmark. The decisive next step is IEEE 39-bus in ANDES with realistic
SG/GFL/GFM/PLL models.

---

## 3. Solid Model Backbone

Base swing / networked second-order model:

```text
M theta_ddot + D theta_dot + L theta = p
```

where:

- `L` is the weighted graph Laplacian, i.e. the linearized active-power Jacobian
  `dP_e / dtheta` under standard lossless assumptions.
- `M = diag(M_i)` is inertia or virtual inertia.
- `D = diag(D_i)` is damping or an equivalent nodal damping approximation.

Inertia-normalized coordinates:

```text
vartheta = M^(1/2) theta
L_tilde = M^(-1/2) L M^(-1/2)
D_tilde = M^(-1/2) D M^(-1/2)

vartheta_ddot + D_tilde vartheta_dot + L_tilde vartheta = M^(-1/2) p
```

Quadratic eigenvalue problem:

```text
(s^2 I + s D_tilde + L_tilde) v = 0
```

If `L_tilde` and `D_tilde` commute, modes decouple. If not, modes couple and
the diagonal formula is only approximate. This is the gap C5 exploits.

Damping margin:

```text
zeta_min = min_k -Re(lambda_k) / |lambda_k|
```

over the complex eigenvalues of:

```text
A = [[0, I],
     [-M^(-1) L, -M^(-1) D]]
```

Important: `nu_2` or `lambda_2(L_tilde)` is a frequency / stiffness measure, not
a damping margin. Higher `nu_2` does not automatically mean safer.

---

## 4. C5: Damping-Selected Reduced QEP

Algorithm:

```text
1. Form L_tilde and D_tilde.
2. Eigendecompose L_tilde -> (nu_k, q_k).
3. Score every mode:
   zeta_hat_k = (q_k^T D_tilde q_k) / (2 sqrt(nu_k)).
4. Select the r modes with the smallest zeta_hat_k.
   Key point: select by damping, not by frequency.
5. Let V contain those selected q_k columns.
6. Build reduced operators:
   L_r = V^T L_tilde V
   D_r = V^T D_tilde V
7. Solve the reduced QEP:
   (s^2 I_r + s D_r + L_r) w = 0
8. Return the smallest damping ratio among the complex roots.
```

Validated so far:

- random systems up to 118 buses,
- `r` around 5-8,
- median error about 0.5-2.5 percent,
- large speedups over full eigensolve,
- frequency-based mode selection fails much more often.

Status:

- Verified in random synthetic settings.
- Originality not fully confirmed. Literature search still required for the
  exact claim: damping-selected reduced QEP for damping-margin estimation under
  non-proportional inverter damping.

---

## 5. C7: Mode-Line Sensitivity Map

For an edge `(i,j)` with weight `w_ij`, the sensitivity of inertial-Laplacian
mode `nu_k` is:

```text
S[k,(i,j)] = d nu_k / d w_ij
           = (q[k,i]/sqrt(M_i) - q[k,j]/sqrt(M_j))^2
```

Interpretation:

- rows: modes,
- columns: lines,
- entries: how much reinforcing a line raises a specific modal frequency.

Validated so far:

- in synthetic tests, the line that most affects the critical mode is often not
  the line that most affects the next mode.

Status:

- The Rayleigh sensitivity formula is classical.
- The mode-resolved planning use may be a contribution, but originality is not
  confirmed until a targeted literature search is completed.

Warning:

- This map raises modal stiffness/frequency, not necessarily damping margin.
  For stability planning, compare against actual or C5-estimated `zeta_min`.

---

## 6. Dead Hypotheses: Do Not Revive Without New Evidence

These were numerically or mechanistically refuted.

| ID | Hypothesis | Why it failed |
| --- | --- | --- |
| H1 | Dramatic transient growth from PLL non-normality | At equal damping ratio, PLL did not add the claimed transient growth. The apparent effect came from low damping. |
| H2 | Weak pathways are generic and robust | Artifact of a toy PLL model. Rare in the corrected 7-state model. |
| H3 | Commutator norm `||[L_tilde,D_tilde]||` is a fragility metric | Did not transfer to the corrected GFL model. |
| H4 | Dynamic Braess paradox from line reinforcement | Cases were explained by reduced modal damping, not a topological paradox. |
| H5 | Clean graph-frequency conversion law | Failed quantitatively; use only as descriptive language if needed. |
| H6 | Resonance-PLL collision `sqrt(K_i) ~= sqrt(nu_c)` predicts instability | Poor stability-boundary prediction and low correlation. |
| H7 | Passivity loss coincides with graph resonances | Non-passive band was mostly near zero frequency, not at graph resonance. |
| H8 | `L_D = D^(-1/2) L D^(-1/2)` is an independent second geometry | Redundant for the same system; GFL also lacks true nodal `D_i` in realistic models. |
| H9 | Commutator norm predicts proportional-formula error | Correlation was near zero. |
| H10 | `nu_c` itself is the stability margin | False. `nu_c` is frequency/stiffness; damping margin is `zeta_min`. |
| H11 | Hidden operational margin / OPF framing | This belongs to steady-state operational planning, not the present dynamic-stability scope. |

Meta-lesson: a high number does not prove a mechanism. Always verify the
mechanism with the right control.

---

## 7. Model-Fidelity Notes

Toy PLL model:

- injected PLL-like terms directly into stiffness,
- not a real PLL,
- use only for intuition, not evidence.

Corrected 7-state 3-node model:

- includes SG, GFM-like dynamics, and GFL PLL states,
- better for mechanism testing,
- still not enough for final claims.

Required next model:

- IEEE 39-bus in ANDES,
- realistic SG/GFL/GFM/PLL models,
- power flow and small-signal linearization,
- compare C5/C7 against full eigenanalysis and time-domain simulations.

ANDES matters because it supports power flow, time-domain simulation, eigenvalue
analysis, symbolic-numeric modeling, and renewable-energy models.

---

## 8. Priority Next Steps

1. Validate on IEEE 39-bus in ANDES.
   - Load case.
   - Run power flow.
   - Build or extract linearized small-signal model.
   - Include realistic SG/GFL/GFM/PLL models.
   - Compare full eigenvalues against C5 reduced-QEP estimates.
   - Test several IBR penetration levels.

2. Validate C7 against full eigen-sensitivities.
   - Do finite-difference line perturbations.
   - Compare predicted modal changes against actual modal changes.
   - Separately evaluate effects on `nu_k` and `zeta_min`.

3. Test single and multi-node SG-to-IBR conversion.
   - Single-node conversion is promising in synthetic tests.
   - Multi-node progressive decarbonization is not yet tested.
   - Check interaction effects.

4. Literature search before novelty claims.
   Suggested queries:
   - "reduced-order damping ratio estimation inverter"
   - "non-proportional damping power network modal reduction"
   - "graph spectral damping margin IBR"
   - "quadratic eigenvalue damping ratio power grid inverter"
   - "mode line sensitivity damping inter-area oscillation"

5. Optional high-risk extensions.
   - Delay-dependent margin using reduced critical-mode operators.
   - Disturbance-shape hidden margin using singular values of the reduced
     resolvent.
   - Sheaf / vector-bundle graph models for heterogeneous SG/GFM/GFL nodes.
   - Hopf/Krein-style mode-collision diagnostics for multi-time-scale IBR
     systems.

---

## 9. Honest Paper Framing

Possible pitch:

> A computationally cheap spectral framework for planning the SG-to-inverter
> transition: estimate damping margin with a damping-selected reduced QEP,
> identify safe conversion and compensation locations, and map which line
> reinforcements target which oscillation modes.

Be precise:

- dispersion laws and Rayleigh sensitivities are classical,
- C5's damping-based mode selection is the candidate main contribution,
- conversion ranking is a strong application,
- C7 is useful but must be framed as mode-resolved planning, not as a new
  eigenvalue formula,
- claims are not final until IEEE 39-bus / ANDES validation is done.

Do not reuse old abstracts that rely on refuted claims such as weak pathways,
large transient amplification from non-normality, or dynamic Braess paradox.

---

## 10. Reporting Template

Every result reported to the user should use this structure:

```text
Claim:

What I verified:

Compared against:

Mechanism checked:

Alternative causes ruled out:

What I did not test:

Originality status:

Confidence:
```

If the result cannot fill this template, it is not ready to be presented as a
paper result.
