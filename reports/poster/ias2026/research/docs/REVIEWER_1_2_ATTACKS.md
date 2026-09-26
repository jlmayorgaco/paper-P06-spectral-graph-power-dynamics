# Reviewer 1 / Reviewer 2 attacks on the TX4 manuscript and evidence-backed responses

- **Manuscript:** `reports/papers/tx4_policy_dependent_incompatibility/main.pdf`
  (11 pp, IEEEtran).
- **Evidence freeze:** tag `IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE` at
  `b9f274e2`.
- **Post-freeze evidence addenda cited by the paper:**
  - `f775db89`: FC18 boundary anatomy;
  - `7a772808` and `f2946257`: Monte Carlo preregistration and amendment;
  - `cfe9fe4b`: Monte Carlo results.
- **Scope of this document:** it contains no new experiments. Every response
  cites existing evidence. The **Concede** lines state what the paper does
  *not* defend.

---

## 1. "This is just generalized Nyquist applied portfolio by portfolio."

**Response.** Yes for the per-portfolio test, and the paper says so (Sec. I,
"What is not claimed"; Table VII, first row). The contributions are
object-level:
- **C1:** the policy-dependent minimal-coalition structure `H(theta)`, its
  exact transverse definition, and the any-order-safety equivalence.
- **C2:** where the non-additivity comes from, namely network closure of local
  actions, and at which interaction order each boundary closes.
- **C3:** exact boundary normals and zero-frequency closure.

GN answers "is this configuration stable". It answers neither "which minimal
coalitions fail" nor "why this coalition and not its subsets".

**Evidence.**
- F10: flagship-only Nyquist, single-port impedance and SCR screens miss the
  coalition changes.
- Per-subset GN reproduces the labels exactly, which is why GN is credited as
  the test.

**Concede.** No new stability test is claimed.

## 2. "Minimal failing sets already exist in N-k analysis."

**Response.** N-k minimal defective sets (Bienstock–Verma 2010) are a *static*
(feasibility or cascading) analogue. `H(theta)` differs in three ways:
- it is dynamic (small-signal, transverse);
- it is indexed by a *continuous control policy*, which reorganizes it at
  fixed network and dispatch (30/36/16 distinct hypergraphs on three planes);
- it comes with an exact algebraic anatomy (Theorem 3: principal minors of the
  closure operator `Q`).

**Evidence.** Table VII; Sec. V-B; Fig. 2.

**Concede.** The combinatorial notions (clutter, antichain, minimal non-face)
are classical and credited.

## 3. "Everything depends on the harmonized IEEEX1 exciter."

**Response.** The harmonization (regulator gain and exciter time constant, with
KE and rate feedback dropped) is declared as a model change (Sec. V-A, Table I).
There are two reasons for it:
- the documented IEEEX1 data are themselves unstable in ANDES (six real modes
  near +1.03 s^-1);
- the documented-data variant is therefore not a usable baseline.

The policy coordinates `k` (AVR gain scale), `t` (time scale) and `h`
(heterogeneity) sweep the excitation over wide ranges, and the structure
persists across them (F7A/B/C).

**Evidence.**
- FC14 (ANDES eigenvalues with the IEEEX1 data);
- `docs/F1B_EXCITER_PROVENANCE.md`;
- F7 maps.

**Concede.** A different exciter model could move the boundaries. The claims
are conditional on the declared model, and the abstract, Sec. V-A and Sec. VII
say so.

## 4. "With governors your flagship result disappears."

**Response.** Correct, and the paper makes this central rather than hiding it
(abstract; Sec. V-B "Primary frequency control"; Table II; Fig. 3). With the
source's own documented TGOV1N governors (untuned):
- the P4 four-unit portfolio becomes stable (alpha_perp +0.127 → −0.0745) but
  remains the least-damped portfolio;
- policy-dependent incompatibility persists: 10 distinct `H`, and 74/399
  non-empty plane points with kappa = 2, 3 and 4.

This *strengthens* the general claim: `H` is a property of system and policy,
not a permanent label on buses.

**Evidence.** FC03 (`FC03_summary.json`, `FC03_points.csv`).

**Concede.** The specific four-unit witness is conditional on the absence of
primary frequency control.

## 5. "No converter current limits, ride-through or DC link: the nonlinear claims are meaningless."

**Response.** The paper makes no nonlinear composability claim beyond scope.
- It declares guards (Appendix D) and reports that 84 of 96 finite thresholds
  are guard activations.
- It states the result as a negative scope check: inside the guards the
  nonlinear hypergraph equals the spectral one (Sec. VII-A).
- The 684–792 MVA support figure is labelled an envelope requirement, not a
  stability requirement.

**Evidence.** FC05, FC06, FC12; Fig. 10; Table V.

**Concede.** Behaviour beyond the guards is undecided. A limiter-equipped model
is future work and must be separately versioned.

## 6. "No independent implementation of the converter model."

**Response.** Correct, and stated (Sec. VII-B).
- ANDES reproduces the network and power flow, the base inter-area mode
  (within 4.1 %) and the governor damping direction (−0.222 → −0.770 s^-1).
- It does not implement the L0 GFL.

Internal validation uses independent *derivative paths* instead:
- direct solves versus the common realization (P3: 3.2e-13);
- direct finite differences versus port derivatives (normal angle ≤ 4.3e-4°);
- the full nonlinear DAE versus eigenvalues (P4: 80/80, r = 0.969).

**Concede.** No portfolio result is claimed as independently reproduced.

## 7. "Your preregistered P2 test failed."

**Response.** Yes. It is reported as **fail as implemented** in Table VI, in
the text, and in the abstract ("one disclosed rule failure").
- **Cause:** the bisection used the three-state classifier label, so it stopped
  at the classifier's unresolved band. There |Re λ| = 4e-6 to 7e-5, above the
  frozen 1e-6.
- **Post-hoc re-bisection:** on the sign of alpha_perp, with the same paths,
  brackets and threshold, it locates 120/120 events at |Re λ| ≤ 2e-9, all as
  complex crossings at 0.41–0.72 Hz. It is reported *separately* and never
  counted as a preregistered pass.
- **What P2 tests:** Theorem 2, which is a continuity argument. The failure is
  a check-implementation defect, not a counterexample to the theorem.

**Evidence.** MC02 P2; MC04; `configs/ias2026/paper_mc_validation_v1.yaml`.

**Concede.** The preregistered rule as written failed.

## 8. "1.5 % of policy-map cells are unresolved: the maps are not certified."

**Response.**
- The classifier has three outcomes and never assigns an unresolved label
  (Sec. III, "Numerical labels").
- The unresolved cells are oscillatory pairs within the Jacobian-assembly error
  near boundaries, not small eigenvalues.
- Re-auditing all 345,229 points in transverse coordinates changed **no**
  resolved label.
- Unresolved cells are drawn separately in Fig. 2.

**Evidence.** FC01 (`TSQ_ieee39_reaudit.csv`): 5,128 unresolved, 0 resolved
changes.

**Concede.** Labels are numerically exact (`d_axis > 10 eps`), not
interval-certified. The paper says so.

## 9. "Transverse stability is not stability; the real system must restore frequency."

**Response.**
- Frozen benchmark models (IEEE-39, Kundur, IEEE-68) have `D = 0` and no
  governors, so their linearizations carry an exact Jordan double zero.
  "Stable" is undefined without a declared quotient.
- Theorem 1 separates the rotation (a gauge) from the drift (the physical
  absence of primary control). It proves the quotient exact for both the
  linearization and the nonlinear reduced field.
- P1 verifies the nonlinear symmetries to 2.4e-14.
- The governed model, where only the rotation is quotiented, is analysed
  explicitly (Item 4).

**Concede.** Absolute frequency restoration is outside the frozen benchmark.
This is stated as a mandatory scope sentence.

## 10. "'Only the full-order term closes the determinant' is either trivial or an artifact of the expansion."

**Response.**
- **Not trivial.** Theorem 3 guarantees only that the zero lies in
  `det(I + Q_SS)`. It does *not* say which order carries it. A triple boundary
  could have been closed essentially by pair terms, with a negligible triple
  coefficient.
- **Measured.** At every one of the ten frozen boundaries, truncation at order
  |S| − 1 leaves 0.29–0.40 (triples) or 0.13–0.16 (four-unit), against a
  closed value of ≤ 4e-13. At the four-unit boundary the orders 2, 3 and 4
  have comparable magnitude (0.51, 0.45, 0.16).
- **Not an artifact.**
  - The Möbius coefficients of `det(I + Q_SS)` are invariant under
    block-diagonal port-basis changes (Theorem 3; S6: invariance 8.2e-12).
  - Individual split minors are *not* invariant (S6 control: median change
    150 %), which is why only block-touching sums are reported.
  - The local factor is factored out explicitly, and the unnormalized
    `det(I + M_SS)` shows the same pattern (Table IV).
- **Wording.** The statement is kept exactly as observed: "no interaction
  truncation below the cardinality of the changing minimal coalition reaches
  the characteristic zero."

**Concede.**
- The expansion is relative to the base (empty) portfolio and to the chosen
  splitting of local and network terms.
- It is an exact bookkeeping of the characteristic function, not a causal
  decomposition.
- For pairs, the full order is order 2, so the pair rows are consistent but
  uninformative.

## 11. "Why not just use displaced MW, SCR, or inertia?"

**Response.**
- **Empirical.** Among 25 four-unit portfolios matched on displaced dispatch,
  all keep the inter-area family damped, while 8 of 10 failing portfolios have
  it unstable (Fig. 4). Displaced dispatch predicts the minimal coalition with
  AUC 0.57. SCR, impedance and flagship-only screens miss coalition changes
  (F10).
- **Structural (Sec. IV-B).** Aggregate dispatch and any per-unit index are
  blind to the off-diagonal network-closure blocks `K_o`. Theorem 3 shows that
  those blocks carry the *entire* non-additive part.

**Evidence.** FC04, UC02, F10.

**Concede.** Rating and inertia can still rank risk. They do not recover the
witness structure.

## 12. "Where are the physical paths or cycles that cause the failure?"

**Response.** Deliberately not claimed, because both hypotheses were tested and
rejected.
1. A bus-diagonal Neumann walk expansion of `T_0^{-1}` does not converge in the
   inter-area band: rho(D_g^{-1} N) = 1.003–1.556 over 0.3–1.5 Hz, both at the
   four-unit boundary and for a Pg-matched stable control. No n-hop
   attribution is therefore claimed (Sec. V-C, "Limitation").
2. The simple-cycle (holonomy) explanation was falsified earlier (E18; the
   dispatch-corrected rerun gives ratio 0.72).
3. The effective rank of `Q` stays 7–8 while kappa changes 4 → 3 → 2 → 3 → 4,
   so interaction dimension is not the mechanism either (Fig. 6, negative
   control).

**Concede.**
- The mechanism is stated at the level of interaction *order*.
- Whether another splitting yields a convergent, interpretable expansion is
  open; the manuscript does not pursue it.

## 13. "IEEE-68 shows nothing, so the phenomenon is an IEEE-39 artifact."

**Response.**
- The IEEE-68 negative is reported, not hidden: the documented model is
  reproduced to 15/15 modes within 5e-4 Hz, and the preregistered
  four-candidate map is empty.
- Policy dependence appears there only at ≥ 6/12 plants
  (kappa 6 → 9 → 11).
- Kundur shows policy dependence on 61/61 and 36/36 lines.
- The theory (C1–C3) is system independent. The *structure* is system
  dependent, and the paper says so ("not universal").

**Concede.** The rich witness structure is demonstrated on IEEE-39 (and
Kundur), not across benchmarks.

## 14. "The lattice is exponential; this does not scale."

**Response.**
- Agreed in general, and stated (Sec. VII-B; complexity statement, credited to
  the sparse-PCA reduction).
- The factorization (7) evaluates every portfolio from one base network and m
  local 2×2 increments at a common operating point. This reduces the *cost per
  label*, not the number of labels.
- For m = 4 the lattice has 16 portfolios; the census has m = 9 (512).
- The monotone class (where pruning is valid) is shown to exclude IEEE-39, so
  no pruning is claimed.

**Concede.** No scalable algorithm for large m is claimed.

## 15. "What is the planning relevance?"

**Response.** It is modest and stated as such (Sec. V-E, Table V).
- Any-order safety is exactly equivalent to hyperedge-free (Proposition 2),
  which gives a set-packing constraint.
- Order dependence exists but is rare (3/327 census targets; 131/30,755 on
  random paths).
- The MW-optimal census plan is any-order safe.
- The policy moves kappa between 1 and infinity at fixed network and dispatch.
- The port derivative gives exact retuning directions, and Newton converges in
  5 steps.

**Concede.**
- Sequencing is not the dominant application.
- A single linear retuning step is not accurate (30 % of the required change).

## 16. (Additional) "Matched dispatch and the common realization are artificial."

**Response.**
- Matched dispatch keeps the AC power flow identical across portfolios, which
  isolates the dynamic effect of replacement from redispatch.
- The common realization (both devices at every candidate, with stable fillers)
  is an exact device for the algebra.
  - Its vertex spectra equal the actual portfolios plus the fillers' −1
    eigenvalues.
  - Direct-path solutions agree with it (P3: 3.2e-13).
  - The direct-path locality residual is ≤ 8e-8.

**Concede.** Redispatch-coupled replacement is not analysed.

## 17. (Additional) "Monte Carlo on your own model is self-validation."

**Response.**
- The synthetic suites test the theorems on systems *unrelated* to the
  benchmark, including expected-failure controls:
  - the magnitude cutoff fails in 235/235 positive planted modes;
  - the non-Metzler control violates monotonicity in 997/1000 families.
- The physical suites test *different code paths*: direct versus common
  realization, finite differences versus port formula, and full nonlinear TDS
  versus eigenvalues.
- The pass rules were committed before any draw, and the one amendment was
  made before the run it governs.

**Concede.** Several synthetic checks are logical identities that validate the
implementation, and the paper labels them so.
