# TX3 Final Reviewer 2 Audit

## Manuscript information

- **Title:** *Finite Connected Externalities in Converter-Rich Power-System Dynamics: Re-Equilibrated Reconstruction, Pole Anatomy, and the Boundary of Modal Materiality*
- **Review date:** 2026-08-28
- **Review round:** Final verification review
- **Role:** Peer Reviewer 2 -- power-systems domain
- **Identity:** Power-system dynamics researcher specializing in converter-rich systems, small-signal stability, differential--algebraic modeling, modal analysis, and operating-point feasibility.

## Review focus

This terminal review evaluates the current PDF's power-system relevance, physical magnitudes, AC re-equilibration, modal-family interpretation, operational materiality, and comparison with conventional sensitivity methods. No additional experiment is considered necessary or requested.

## Overall assessment

### Recommendation

- [x] **Accept**
- [ ] Minor Revision
- [ ] Major Revision
- [ ] Reject

### Confidence

**4/5 -- High confidence.**

### Summary assessment

The manuscript develops a finite connected intervention calculus in which every subset of a controller, dispatch, or network-action portfolio receives a new AC solution, dynamic initialization, and descriptor model. It demonstrates direct-versus-connected reconstruction, exact Schur-factor anatomy, reproducible finite-pole externalities, and a separate test of whether those effects reach the actual limiting modal family.

The current version is materially stronger from a power-systems perspective. It identifies the bus-39 slack policy, fixed-P/Q GFL representation, system and device bases, controller and network endpoints, voltage/current ranges, and implemented limits. It also reports the leading synchronous-machine states of the 1.22-Hz critical mode and distinguishes the finite-coalition estimand from first-order sensitivity, higher-order local sensitivity, impedance analysis, and proximal eigenvalue displacement. Most importantly, it preserves the negative findings: nominal and load-conditioned materiality are rejected, while E05D remains unresolved because all eight empty baselines were ineligible and no coalition result was inspected. Remaining concerns concern wording, provenance, and presentation rather than scientific validity. The manuscript does not convert externality existence into a stability hazard or industrial mitigation claim. I recommend acceptance without further experimental work.

## Strengths

### S1: The AC re-equilibration policy is electrically explicit

Section III-A, p. 3, identifies the 100-MVA AC base, bus 39 as the reference/slack bus, GFL buses 36--38 as fixed-P/Q injections, disabled legacy PV voltage equations, retained synchronous PV controls, and slack balancing of A7. It also reports a 0.84056-pu minimum voltage, 0.90842-pu maximum GFL current, and a 1.2-pu current limit. This resolves the principal ambiguity about how finite action vertices are physically closed.

### S2: Action magnitudes are interpretable in engineering units

Table II and Section III-B, p. 4, report A7 as 540 to 560 MW, A8 as 0.0232 to 0.02552 pu, 1000-MVA converter ratings, baseline GFL injections of 560/540/830 MW, PLL frequency 2.0 to 2.4 Hz, and current-command time constant 0.0200 to 0.0167 s. Readers can assess the scale of the frozen actions rather than seeing only percentages.

### S3: The novelty perimeter is clear and restrained

Table I, p. 3, distinguishes first-order eigen-sensitivity, higher-order total sensitivity, return-ratio/impedance analysis, proximal eigenvalue displacement, M\"obius coefficients, and the present finite re-equilibrated workflow. Sections I--II explicitly state that the contribution is not a new M\"obius identity, Schur complement, or stability criterion.

### S4: Structural conclusions are supported by independent closures

Sections IV--V, pp. 5--6, retain 266 resolved pair cases, 147 resolved triple cases, 129/129 locked-holdout reconstructions within five uncertainty budgets, and Schur-factor residuals no larger than \(7.39\times10^{-13}\). The 6.3-fold A7--A8 frozen-affine underestimate is limited to OP00 rather than promoted as a universal correction.

### S5: Modal consequence and negative evidence are handled correctly

Section VI, pp. 7--8, preserves the rejected nominal and load-conditioned materiality findings and the **unresolved--baseline ineligible** E05D status. The critical-mode audit reports a synchronous-electromechanical family, leading GENROU32/30/34/31/33 angle channels, 200/200 successful subset tracks, and a maximum critical-mode damping externality of only \(5.14674\times10^{-6}\). This supports a family mismatch without creating a post hoc gate or intervention claim.

## Weaknesses

### W1: Slack-policy dependence and the two voltage contexts could be stated more explicitly

**Problem:** Section III-A reports a minimum E03 vertex voltage of 0.84056 pu, while Section III-D separately defines the E05C endpoint by the first 0.85-pu voltage boundary. The PDF does not explicitly distinguish the campaign populations or note that the algebraic A7--A8 result is conditional on centralized bus-39 balancing.

**Why it matters:** Readers may infer an inconsistency or treat the algebraic externality as independent of balancing policy.

**Suggestion:** Add one sentence distinguishing E03 and E05C voltage populations and acknowledging the registered slack allocation.

**Severity:** Minor.

### W2: “Coherent multi-machine mode” is stronger than the displayed evidence

**Problem:** Section VI-D calls the limiting family a “coherent multi-machine electromechanical mode,” while the PDF reports unsigned descriptor-endogeneity/participation shares but no relative rotor-angle phase information.

**Why it matters:** Coherency normally implies a phase relationship among machine motions, not merely distributed participation.

**Suggestion:** Use “distributed multi-machine electromechanical mode,” or report the retained phase relationship.

**Severity:** Minor.

### W3: Benchmark and model provenance remain lightly cited

**Problem:** Section III-A cites ANDES, but the reference list does not identify the provenance of the precise IEEE 39-bus realization, conversion of buses 36--38 to GFL plants, or detailed REGCP/PLL model definitions.

**Why it matters:** The package supports reproducibility, but readers should identify the model lineage directly from the article.

**Suggestion:** Add citations or a precise package pointer for the system realization and converter-model definitions.

**Severity:** Minor.

## Detailed comments

### Title and abstract

The title accurately describes the finite, re-equilibrated scope and materiality boundary. The abstract uses “registered damping margin,” avoiding a universal standard, balances structural successes with the negative modal result, and explicitly rejects stability or mitigation certification.

### Introduction and literature

The research question is clear. Table I improves comparison with established methods, and the manuscript correctly states that impedance, return-ratio, modal sensitivity, and M\"obius inversion answer related but non-identical questions. Benchmark and GFL-model provenance are the only notable domain-literature gaps.

### Methodology

The DAE, action-hypercube, uncertainty, holdout, factorization, contour, and mode-tracking protocols are well specified. The fixed-P/Q/slack description makes re-equilibration reproducible in physical terms. Action endpoints and ratings are concrete. The materiality thresholds remain correctly identified as preregistered study criteria rather than grid-code limits.

### Results

Direct and connected calculations have distinct numerical paths and close case by case. The A7--A8 algebraic-dominant and A3--A6 pole-dominant examples prevent a single-mechanism interpretation. A3--A6's failed contour certificate remains visible. E05D is interpreted correctly: baseline ineligibility supports neither acceptance nor rejection of weak-grid materiality. The critical-mode audit remains descriptive.

### Discussion and conclusion

The revised sensitivity subsection accurately refers to a baseline first-order screen. The paper clearly states that a dimensionless determinant externality is not damping, instability probability, energy, or economic cost. Limitations exclude industrial mitigation, cycle/SCC causality, nonlinear consequence, and blind EMT transfer. The conclusion preserves C1, C2, D1, D2, conditional D3, the nominal/load rejections, and weak-grid unresolved status.

### References

Citation formatting is consistent and readable. The main theoretical comparators are covered. Benchmark-specific and model-provenance citations would improve traceability but do not affect acceptance.

## Questions for the authors

1. Can the authors clarify why the E03 minimum voltage of 0.84056 pu and the E05C 0.85-pu endpoint refer to distinct registered stages?
2. Is “coherent” based on retained rotor-angle phase relationships, or is “distributed multi-machine” intended?
3. Can the exact IEEE 39-bus/GFL realization be named directly in the manuscript or package pointer?

## Minor issues

- Capitalize “AC” consistently.
- Consider replacing “coherent” with “distributed” unless relative machine phases are reported.
- Table I is informative but uses small type.
- Fig. 3's triple terms are visually difficult to distinguish from zero; the numerical text compensates.
- Table IV contains tightly wrapped cells.
- The 10-page PDF is clean, balanced, and free of visible clipping or overlapping elements.

## Dimension scores

| Dimension | Score | Descriptor | Notes |
|---|---:|---|---|
| Originality (20%) | 84 | Strong | Clear finite re-equilibrated portfolio estimand and integrated evidence ladder. |
| Methodological rigor (25%) | 88 | Strong | Registered protocols, independent closures, explicit AC policy, and reproducible tracking. |
| Evidence sufficiency (25%) | 83 | Strong | Extensive internal validation and negative evidence; external validity is intentionally narrow. |
| Argument coherence (15%) | 93 | Exceptional | Clear progression from existence through anatomy to materiality boundary. |
| Writing quality (15%) | 91 | Exceptional | Precise, conservative, and well structured. |
| Literature integration | 78 | Strong | Strong comparator coverage; benchmark/model provenance could be expanded. |
| Significance and impact | 79 | Strong | Important diagnostic distinction with bounded practical reach. |
| **Weighted average** | **87.2/100** | **Accept** | No scientific defect requiring further experimentation or re-review. |
