## 2. Research questions

The campaign tests one meta-hypothesis and asks whether it survives:

> static weakness ≠ dynamic weakness ≠ contextual weakness ≠ intervention leverage.

The question is not "which bus is weakest". It is which properties of nodes,
links, corridors and control locations robustly determine, or move, the dynamic
incompatibility boundaries when any of the following change:
- the replacement portfolio;
- the controller policy;
- the operating point;
- uncertain parameters.

The preregistered hypotheses (`docs/CDW_PREREG_V1.md`) are:

| id | hypothesis |
|---|---|
| H1 | contextual sign reversal recurs |
| H1t | the reversal survives mode tracking |
| H1b | fixed node-only rankings are insufficient |
| H1c | contextuality is robust under envelopes |
| HS | submodularity or supermodularity |
| H2 | frozen vs total sensitivity; dynamic vs static baselines |
| H3 | robust weak corridors |
| H4 | topology alone |
| H5 | plan-level design |
| H6 | spectral closure and modal mixing |
| H7 | certified reduction |
| H8 | modal energy vs the limiting mode |
| H9 | cross-model transfer |

## 3. Relation to the frozen TX4 work

- **Frozen and read only.** TX4 is frozen at tag `TX4_FINAL_MANUSCRIPT_FREEZE`
  (69f200df). The CDW branch starts from that tag, and no TX4 file, result,
  ledger or tag was modified.
- **Inherited theorems** (CDW theory §A):
  - the transverse quotient;
  - the any-order-safety characterization of the minimal incompatibility
    hypergraph;
  - boundary localization;
  - exact network-closure factorization;
  - the boundary-sensitivity formula;
  - the deflated zero-frequency port.
- **Inherited IEEE-39 evidence (CDW theory §B).** It is benchmark-specific and
  is used as motivation only:
  - the nominal P4 witness;
  - the policy dependence of H;
  - the one-condition line ranking;
  - the robustness of the phenomenon, but not of the witness identity.
- **The TX4 ParaEMT line is closed.** No CDW statement is an EMT result.

## 4. Graph-DAE model

The model is the frozen TX4 IEEE-39 phasor DAE (`theory/CDW_GRAPH_DAE_V1.md`):
- ten two-axis SGs (bus 39 is the interconnection);
- the frozen GFL, with SRF-PLL, PI outer and inner loops, and a leaky Q/V
  regulator whose gain is `g`;
- constant-power loads;
- an exact branch-terminal network, with 12 off-nominal-tap transformers, line
  charging and 2 shunts.

**Network.** The admittance is `Y(ρ) = Y_sh + Σ_e ρ_e C_eᵀ Y_e C_e`, where `ρ_e`
scales the whole two-port. A unit test checks that `Y(1)` equals the frozen
admittance matrix bit for bit.

**Graph operator.** The coupling Laplacian `L_B` has weights `b_e/t_e`, and the
decomposition `Im Y = −L_B − Δ_tap + B_ch + B_sh` is exact (unit test). No
`B diag(y) Bᵀ` simplification is used.

**Replacement.** Replacing SG `i` by a GFL keeps matched dispatch and an equal
rating, so the AC operating point is identical for every portfolio, as in TX4.

**Policy.** `θ = (g, k, t, h)` collects the Q/V gain and the scale, time-constant
scale and heterogeneity of the AVRs.

## 5. Definitions of weakness

Five distinct objects are defined (theory §C3). Their coincidence is tested, not
assumed.
- **W_static.** SCR, Thevenin |Z|, gSCR, electrical distance, effective
  resistance, Fiedler entries and edge scores, betweenness, flows, dV/dQ.
- **W_dynamic.** The boundary sensitivity `dRe s*/da`, with its sign, in both a
  frozen-operating-point and a re-equilibrated total version.
- **W_contextual.** The distribution of `Δ_i α(S) = α(S ∪ {i}) − α(S)` over
  contexts `S`. It is a property of an (intervention, context) pair, not of a
  bus.
- **W_closure.** The split of the boundary derivative through the exact closure
  factor into device, driving-point and transfer-coupling parts. Physical
  ablations are allowed; zeroing blocks of `K_o` is not used.
- **W_control.** The finite leverage within the admissible ranges.

**Frozen vs total** (theory §C4). Let the operating point be `w = (x, z, r)`:
states, voltages and device references. Then
`dw*/da = −R_w^{-1} R_a` and `dλ/da = v^H (∂_a A + D_w A[dw*/da]) u / v^H u`.

Two equilibrium semantics are used:
- **SPR** (schedule-preserving re-initialization, the TX4 semantics): schedules
  are held and references are re-derived;
- **RP** (reference-preserving): references are held and only the slack
  mechanical power is free.

**Structural remarks.** These are exact consequences of the equations and are
verified numerically in E3.
- Under SPR, the frozen and total derivatives coincide for every
  controller-only coordinate (Q/V gain, PLL, AVR gain).
- Setpoints have a frozen partial that is identically 0.

## 6. Contextual sign-reversal theory

**Definitions.**
- A **contextual sign reversal** of intervention `i` needs two contexts: one
  with `Δ_i α(S1) ≤ −τ` and one with `Δ_i α(S2) ≥ τ`.
  - The material threshold is `τ = 0.01 s⁻¹`, justified in prereg §1.
  - A **stable-context reversal** requires both contexts to be transversely
    stable, which is the planner's situation.
- `α` is a maximum over modes, so a sign change of `Δ_i α` can be produced by a
  switch of the rightmost mode. A **mode-tracked** reversal therefore also
  requires two things:
  - the rightmost mode of each context must be matched, by bus-voltage mode
    shape MAC ≥ 0.8 and |Δf| ≤ 0.15 Hz, to a mode of `S ∪ {i}`;
  - that matched mode must be the rightmost mode of `S ∪ {i}` (same-mode
    transitions).

**Structural facts.**
- On the lattice `2^V`, submodularity of a set function is equivalent to the
  one-step condition `α(S∪i∪j) − α(S∪i) − α(S∪j) + α(S) ≤ 0` (unit test). This
  makes every violation a minimal counterexample of cardinality one.
- A fixed node-only ranking is judged by the most favourable ranking possible:
  the exact linear-ordering optimum over the policy's own resolved context
  triples, computed by dynamic programming over subsets. Its failure is
  therefore not an artefact of a poorly chosen ranking.
