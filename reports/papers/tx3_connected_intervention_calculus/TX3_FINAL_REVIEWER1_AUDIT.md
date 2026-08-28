# TX3 Final Reviewer 1 Audit

## Manuscript information

- **Title:** *Finite Connected Externalities in Converter-Rich Power-System Dynamics: Re-Equilibrated Reconstruction, Pole Anatomy, and the Boundary of Modal Materiality*
- **Manuscript ID:** TX3
- **Review date:** 2026-08-28
- **Review round:** Terminal re-review
- **Role:** Peer Reviewer 1 -- novelty and theoretical positioning
- **Identity:** Power-system dynamics researcher specializing in converter interactions, multivariable feedback, impedance-based stability, descriptor systems, perturbation determinants, and higher-order modal sensitivity.

## Review focus

This review evaluates the compiled manuscript against return-ratio and impedance analysis, perturbation-determinant and spectral-shift theory, higher-order sensitivity, M\"obius inversion, and the Dewangan--Puricelli--Beerten proximal eigenvalue-displacement method. The experimental scope is treated as closed.

## Overall assessment

### Recommendation

- [ ] Accept
- [x] **Minor Revision**
- [ ] Major Revision
- [ ] Reject

### Confidence

**4/5 -- High confidence.**

### Summary assessment

The manuscript defines finite pair and triple M\"obius coefficients of a descriptor log-determinant after every required action subset is independently AC re-equilibrated, initialized, and linearized. It reconstructs these coefficients from connected total derivatives, separates algebraic, mass, and finite-pole factors, and then tests whether the identified interactions affect the actual damping-limiting mode.

The novelty is now presented credibly as an integrated, experimentally gated workflow rather than a new M\"obius identity, determinant formula, or sensitivity theorem. Table I gives an explicit method-to-method perimeter, Section V-C distinguishes the estimand from proximal eigenvalue displacement, and the relative-determinant paragraph correctly limits contour claims to finite-dimensional determinant identities and the argument principle. The manuscript also confines the 6.3-fold discrepancy to the OP00 frozen-affine diagnostic and narrows its conclusion to cases where finite portfolio composition is the decision question.

The remaining concerns are minor: Table I compares somewhat heterogeneous regimes under the heading “Amplitude,” the generic versus log-determinant-specific parts of the calculus could be separated more explicitly, and the set-function “dividend” terminology would benefit from its established disciplinary lineage. None requires a new experiment or changes the frozen claims.

## Strengths

### S1: The novelty perimeter is explicit and appropriately bounded

Pages 1--3 state that the contribution is not a new M\"obius identity, Schur complement, trace formula, perturbation determinant, or eigenvalue sensitivity. Table I compares first- and higher-order sensitivity, return-ratio/impedance, proximal eigenvalue displacement, M\"obius coefficients, and the proposed workflow using common dimensions. This materially supports the integrated-workflow novelty claim.

### S2: The closest recent comparator is distinguished at the object level

Section V-C, p. 6, identifies Dewangan--Puricelli--Beerten as an open-loop proximity and first-/second-order eigenvalue-displacement method for identifying converter interaction modes. It then distinguishes the present finite action-coalition estimand, subset-wise AC re-equilibration, descriptor-factor anatomy, and independently selected minimum-damping family.

### S3: The perturbation-determinant boundary is technically responsible

Page 2 explicitly defines the relative determinant

\[
D(s,a)=\det T(s,a)/\det T(s,0)
\]

as a finite-dimensional non-normal determinant ratio. It states that contour results follow from the argument principle and registered numerical certificates, not self-adjoint spectral-shift theory. Pages 6--7 retain the term “connected pole localization,” avoiding spectral-causality overclaim.

### S4: The sensitivity comparison matches the implemented diagnostic

Page 5 confines the 6.3-fold discrepancy to A7--A8 at OP00 and acknowledges that higher-order equilibrium sensitivities can be computed. Section VII-B, p. 9, is titled “Finite Portfolio Information Absent from a Baseline First-Order Screen,” accurately identifying the comparator used.

### S5: The final interpretation remains consistent with the frozen evidence

Pages 7--10 preserve the separation between reproducible dynamic externality and registered modal materiality. Nominal and load-conditioned materiality remain rejected, weak-grid materiality remains unresolved because of baseline ineligibility, and no cycle/SCC, mitigation, or EMT-transfer claim is introduced.

## Weaknesses

### W1: Table I places heterogeneous regimes under “Amplitude”

**Problem:** “Local,” “Local Taylor series,” “Frequency-domain model,” “First-/second-order modal displacement,” and “Finite” appear in one column. A frequency-domain model is not itself an amplitude regime.

**Why it matters:** The table is central to the novelty argument, and a heterogeneous column weakens its precision.

**Suggestion:** Rename the column “Action/perturbation regime” and state that rows summarize representative formulations, not strict limits on all extensions.

**Severity:** Minor.

### W2: Functional-generic and log-determinant-specific contributions are not fully separated

**Problem:** Equations (3)--(4) apply to any sufficiently smooth scalar set function, while the connected trace expansion, Schur anatomy, and contour localization depend on the descriptor log-determinant.

**Why it matters:** “Finite connected intervention calculus” could otherwise be read either too broadly or too narrowly.

**Suggestion:** State explicitly which parts are functional-generic and which require the descriptor log-determinant.

**Severity:** Minor.

### W3: “Set-function dividend” lacks its wider established terminology

**Problem:** Page 2 calls the triple coefficient a “third-order set-function dividend (M\"obius coefficient),” but related Harsanyi-dividend and interaction-index terminology is not acknowledged.

**Why it matters:** The paper correctly does not claim M\"obius inversion as new; a terminology sentence would further prevent “externality” from appearing to rename the combinatorial object.

**Suggestion:** Add the established synonym and a foundational reference, reserving “externality” for the engineering interpretation.

**Severity:** Minor.

### W4: A few compressed phrases remain opaque

**Problem:** The abstract says that cases exceeded “five conservative uncertainty budgets,” and Fig. 1 says cases close “within five budgets.” Table I also uses very small type.

**Why it matters:** Prominent claims should be self-contained before the protocol is read.

**Suggestion:** Use “five times the case-specific numerical uncertainty budget” on first occurrence and shorten Table I entries if possible.

**Severity:** Minor.

## Detailed comments

### Title and abstract

The title is accurate and does not promise unsupported causality or mitigation. The abstract includes the structural result, factor closure, OP00 qualification, critical-mode mismatch, and negative materiality conclusion. “Registered damping margin” is appropriately narrower than a universal operational-margin claim.

### Introduction and literature

The Introduction now establishes the contribution through both positive definition and explicit non-claims. Table I is effective in showing that novelty lies in the registered combination of finite re-equilibrated vertices, connected reconstruction, descriptor-factor anatomy, and a separate modal-materiality test. The treatment of modal sensitivity, return-ratio/impedance, M\"obius inversion, determinant theory, and proximal eigenvalue displacement is sufficient for the bounded novelty claim.

### Methodology

The estimand is unambiguous: all vertices required for each registered pair and triple are evaluated, not every order of the complete \(2^8\) power set. The manuscript correctly distinguishes a finite M\"obius coefficient from a causal multi-device mechanism. The finite relative determinant and argument-principle boundary materially improve theoretical precision.

### Results

Direct-versus-connected reconstruction and Schur-factor closure are strong independent numerical checks. The paper correctly notes that both routes share the same model and action definitions, so closure validates the calculus and implementation rather than physical causality. The 6.3-fold result is adequately bounded, and the unresolved A3--A6 contour is retained.

### Discussion and conclusion

The sensitivity subsection is appropriately narrow and acknowledges higher-order and equilibrium-aware sensitivities. The three-way distinction among finite-action, re-equilibration, and interaction-order error is useful. The conclusion is aligned with C1, C2, D1, D2, conditional D3, and negative findings N1--N4.

### References

The bibliography is compact but adequate for the power-system novelty perimeter. A set-function dividend/Harsanyi reference would complete the terminology boundary.

## Questions for the authors

1. Is Table I intended to describe representative use of each method or a strict capability boundary?
2. Which contributions survive if \(\Phi\) is replaced by another smooth scalar functional, and which require the descriptor log-determinant?
3. Is “set-function dividend” intended as the established Harsanyi/M\"obius dividend, with “externality” reserved for its engineering interpretation?

## Minor issues

- Expand “five conservative uncertainty budgets” at first use.
- Consider “Action/perturbation regime” instead of “Amplitude” in Table I.
- Table I is informative but visually dense.
- Figure 4 is legible, although its rotated coalition labels are crowded.
- Table IV appropriately replaces “negligible” with “below registered gate on limiting family.”

## Dimension scores

| Dimension | Score | Descriptor | Notes |
|---|---:|---|---|
| Originality (20%) | 76 | Strong | Novel integration of established components is differentiated explicitly. |
| Methodological rigor (25%) | 89 | Strong | Registered gates, locked holdouts, independent routes, and retained negatives are exemplary. |
| Evidence sufficiency (25%) | 84 | Strong | Evidence supports all bounded benchmark claims and the materiality boundary. |
| Argument coherence (15%) | 90 | Exceptional | Finite action to factor anatomy to materiality is internally consistent. |
| Writing quality (15%) | 89 | Strong | Precise and conservative, with minor density issues. |
| Literature integration | 72 | Adequate | Closest power-system comparators are well positioned; set-function terminology remains compact. |
| Significance and impact | 76 | Strong | Valuable for finite portfolio studies even though materiality was not established here. |
| **Weighted average** | **85.30/100** | **Accept range; Minor Revision recommended** | Remaining issues concern precision and presentation, not evidence or scope. |
