# Skeptical novelty review: network frequency moments and PLL compensation

Review date: 2026-10-04. Independent review under `ieee-access-academic-quality`. This is a bounded primary-literature check and an inspection of the current experiment artifacts, not an exhaustive priority search. No simulations or model/design changes were performed. The design run was still producing results when inspected; conclusions below deliberately do not assume its eventual outcome.

## Verdict

The plausible contribution is a **model-specific structural identity and its demonstrated design consequence**: finite PLL retuning has a hierarchy of exact low-frequency effects, and its first nonzero effect is aligned with the network synchronization response. The exact scalar constraint can then distribute delay compensation among sites while a separate calculation restores modal damping.

This overlaps substantially with classical type-II tracking, transfer-function moment matching, simple-pole resolvent theory, and network synchronization analysis. Neither the Laurent expansion, the rank-one residue, signed response-area constraints, nor the constrained QP is new machinery. Combining those ingredients is insufficient by itself. I did not locate a primary paper stating this exact heterogeneous SG/GFL, lossy-network coefficient, but the bounded search does not establish priority.

The strongest defensible story is that **preserving low-frequency frequency-response moments and preserving stability are different design tasks, and the network supplies an exact, inexpensive constraint for doing both**. The claim should be conditional on successful independent full-spectrum validation and a clear reason for preserving the selected moment.

## What the actual artifacts support

I inspected `PROTOCOL.md`, the multimode addendum, implementation-repair record, `moments.py`, the multimode corrector, structural audit, and current moment/design tables.

- The implementation retains the exported 174-state hidden system and ten PLL ports. It restores the two rotational identities explicitly before computing the germ. The audit correctly calls this a numerically validated symmetry-restored model, not an interval certificate for independent binary64 matrix entries.
- The bordered recurrence and small-frequency resolvent comparisons support the proposed leading coefficients for all 39 outputs and the selected three input directions. These are correlated outputs of one model, not 117 independent validations or network instances.
- The computed site weights have mixed signs. They are synchronization-residue coefficients, not positive inertia shares, centralities, or unconditional site rankings.
- The one-mode corrector failed and was replaced by a multimode method. The addendum preserves that failure and labels the revised comparison exploratory. This is scientifically appropriate. A finite pole catalog remains incomplete evidence of stability.
- The current comparison contains unchanged gains, sitewise cancellation, and the constrained collective correction. It does not yet contain the most informative ablation: the identical multimode optimizer without the moment equality.

## Closest verified primary sources

| Source and verified access | What already exists | Fair distinction of the present proposal |
|---|---|---|
| F. Paganini and E. Mallada, **Global Analysis of Synchronization Performance for Power Systems: Bridging the Theory-Practice Gap**, IEEE TAC 65(7), 2020. [Published full text](https://mallada.ece.jhu.edu/pubs/2020-TAC-PM.pdf), especially Section VI-A and Lemma 7. | Network synchronization modes, aggregate versus residual bus-frequency dynamics, and perturbations that preserve DC frequency response. | The proposed higher-order finite PLL identities retain the particular full hidden model and heterogeneous exact PLL delays; they do not require a proportional swing/turbine surrogate. Do not imply that DC invariance or isolating the zero mode is new. |
| J. Wilson, A. Nelson, B. Farhang-Boroujeny, **Parameter Derivation of Type-2 Discrete-Time Phase-Locked Loops Containing Feedback Delays**, IEEE TCAS-II 56(12), 886–890, 2009. [Author full text](https://span.ece.utah.edu/uploads/Parameter_derivatrion_DT_PLL_with_delays.pdf), Sections II–V. | Delay-dependent PLL gain selection, with transient pole targets and separate stability/dominance verification. | This is a discrete-time individual-loop construction; the proposed scalar equality distributes a network moment requirement across multiple PLLs. Compare purpose and assumptions, rather than claiming that gain compensation for PLL delay is new. |
| W.-J. Beyn, **An Integral Method for Solving Nonlinear Eigenvalue Problems**, Linear Algebra and its Applications 436(10), 3839–3863, 2012, DOI 10.1016/j.laa.2011.03.030. [Primary manuscript](https://arxiv.org/pdf/1003.1580), Theorem 2.4. | The simple-eigenvalue Keldysh resolvent residue and its left/right-vector derivative normalization. | The power-system identification of the matrices and the resulting parameter-order hierarchy are the application; the singularity theorem is standard. |
| F. Aalipour and T. Das, **Shaping the Transient Response of Nonlinear Systems to Satisfy a Class of Integral Constraints**, Advances in Control Applications 4(3), e110, 2022, DOI 10.1002/adc2.110. [Primary full manuscript](https://arxiv.org/pdf/2012.12493), Section 2. | Signed transient-area constraints, compensating excursions, and controller structures that enforce them while requiring separate stability analysis. | Its energy-storage/load-following problem differs from PLL network response moments. It nevertheless precludes marketing signed-area compensation as a new control principle. |
| B. Salih and T. Das, **Transient Response of Linear Systems Under Integral Constraints**, ASME Journal of Dynamic Systems, Measurement, and Control 142(12), 2020, DOI 10.1115/1.4048106. [Primary author bibliography](https://mae.ucf.edu/TDas/wp-content/uploads/2023/12/TuhinDas_CV.pdf); also cited and described in the Aalipour–Das manuscript above. | Transfer-function structure associated with integral transient constraints and compensator design. | Bibliography verified through the author; the publisher full text was inaccessible in this review. Do not make fine-grained absence claims about its theorems without obtaining it. |
| U. Zulfiqar, V. Sreeram, X. Du, **Finite-Frequency Power System Reduction**, IJEPES 113, 35–44, 2019. [Publisher article](https://doi.org/10.1016/j.ijepes.2019.05.022). | Low-frequency moment matching and modal preservation in power-system reduction. | The proposed identity constrains interventions on the retained model, rather than fitting a reduced model. Moment preservation and modal requirements remain established neighboring ideas. |

These sources are method/theory comparators, not identical numerical competitors. No source inspected here establishes the precise cross-site formula in the proposed SG/GFL model. Conversely, absence of that exact formula does not establish a substantial engineering advance.

## Exact scope and wording

For fixed integral gains, state the result as

\[
H_b(s)-H_a(s)=-s^2\kappa H_0+O(s^3),\qquad
\kappa=\sum_iw_i\left(\frac{\Delta\tau_i}{K_{i,i}}-
\frac{\Delta K_{p,i}}{K_{i,i}^2}\right).
\]

The coefficient is exact for finite parameter changes within the model contract. The displayed expression without the remainder is **not** an exact transfer-function identity at arbitrary frequency. Finite parameter changes also do not make the linearized theorem exact for finite-amplitude nonlinear events.

For a common step disturbance, stable responses satisfy

\[
\int_0^\infty(f_b-f_a)\,dt=0,\qquad
\int_0^\infty t(f_b-f_a)\,dt=\kappa f_\infty.
\]

Changing integral gains instead generally changes the signed area, with the coefficient built from inverse-gain differences. Call that the zeroth **step-transient** moment, equivalently the first impulse-response moment. Normalize by the final response only when it is nonzero and adequately resolved; for disturbances with zero final response the leading formula vanishes and higher orders are needed.

The theorem needs a gain-independent equilibrium/hidden return, nonzero integral gains, an analytic hidden block near zero, a simple rotational zero, nonzero residue denominator, and the stated frequency-output map. Time-integral interpretations additionally need stability and convergent moments. Common finite-window frequency measurement is compatible with the leading difference when lower-order coefficients agree; specify that measurement contract.

The sitewise law is immediate from preserving the scalar low-frequency combination `Kp - Ki*tau`. The less obvious and potentially useful observation is that only one **network-weighted** equality is required, leaving directions available for damping correction. In the present case the signed weights make that freedom different from uniform averaging. Their signs and magnitudes must not be attributed solely to line losses without a controlled comparison.

## Plausible publishable minimum

1. **A complete structural theorem.** Prove the finite-change coefficient and rank-one alignment for the stated class of models, including input/output factors and the gauge assumptions. Present it as a specialization of resolvent/moment theory. Verify that the model class, rather than one numerically restored matrix, supplies the rotational identities.
2. **A useful requirement.** Explain why preserving this signed moment matters for the intended task. It cannot by itself bound nadir, RoCoF, absolute error, settling time, or energy. If no operational use emerges, position it as a limitation of moment-based response matching: identical moments can conceal materially different stability.
3. **The right ablations.** Compare unchanged gains, sitewise compensation, minimum-cost moment-only compensation, the same multimode optimizer without the equality, and the joint method. Use identical gain bounds, pole margin, trust-region budget, and validation. The unconstrained optimizer measures the price of moment preservation; moment-only versus joint isolates the need for stability correction. A stronger optimizer beating a sitewise algebraic rule is not enough.
4. **Independent evidence.** Complete the exact-DDE spectrum check for the headline candidates. Verify linear moment predictions independently, with controlled integration tails if time-domain integrals are reported. Use nonlinear events to assess practical behavior, without claiming exact nonlinear moment conservation. Preserve the failed one-mode and root-collision attempts.
5. **A small transfer test.** Freeze the revised method before testing at another operating point or topology and at heterogeneous delay patterns. Show whether the signed weights and compensation choice change. This supports a reusable design law; one fixed IEEE39 trim supports a narrower case study.

The present material can become a credible engineering contribution if these gates establish that the exact network constraint is useful and compatible with validated damping restoration. Without them, it is an instructive application of classical moment matching and synchronization residues, with a local optimization example. It does not establish a maximum replacement limit, a necessary SG-support floor, globally optimal tuning, or a new general graph-theoretic stability law.

## Suggested central claim

“For the specified SG/GFL network model, finite PLL gain and delay changes obey exact low-frequency moment identities. These identities expose response features that retuning cannot change and provide a scalar constraint for distributing delay compensation across sites while separately restoring modal damping.”

Use the final clause only after independent spectrum validation succeeds. The empirical superiority or generalization claim remains open.
