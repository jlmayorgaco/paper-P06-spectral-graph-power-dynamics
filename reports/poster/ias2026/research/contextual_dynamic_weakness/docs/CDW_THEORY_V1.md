# CDW theory V1 — Contextual Dynamic Weakness in Inverter-Rich Power Networks

Date: 2026-09-12. Branch `research/contextual-dynamic-weakness`, created from tag
`TX4_FINAL_MANUSCRIPT_FREEZE` (commit 69f200df). Not pushed.

**Scope.** This document separates:
- **A.** Inherited TX4 results that are proved.
- **B.** Inherited TX4 evidence that is empirical and benchmark-specific.
- **C.** New definitions.
- **D.** New hypotheses.
- **E.** What is not yet tested.

Nothing in D or E is a result. No IEEE-39 observation is presented as a
theorem. Model notation is in `theory/CDW_GRAPH_DAE_V1.md`.

---

## A. Inherited, proved (TX4; statements as frozen in `reports/papers/tx4_policy_dependent_incompatibility/main.tex`)

### A1. Nonlinear power-system DAE

The system is

    x' = f(x, z; δ, θ, ρ, π),     0 = g(x, z; δ, θ, ρ, π),

with:
- `x` the device states;
- `z` the rectangular bus voltages;
- `δ ∈ {0,1}^m` the replacement (composition) indicators;
- `θ` the controller parameters;
- `ρ` the network (branch) variables;
- `π` the operating-point variables (dispatch, loads, setpoints).

If `g_z` is nonsingular at an equilibrium, then

    A = f_x − f_z g_z^{-1} g_x.

### A2. Exact transverse quotient (TX4 Thm 1, main.tex 282–296; claim V01)

**Hypotheses.** Rotation covariance holds, together with `D = 0`, constant
`P_m` and a frequency-independent network.

**Statement.**
- `A R_x = 0` and `A w = ω_B R_x`.
- `C = span{R_x, w}` is A-invariant.
- `det(sI − A) = s² det(sI − A_⊥)` with `A_⊥ = Zᵀ A Z`, for any orthonormal
  `Z` whose range is `C^⊥`.
- `α_⊥ = max Re σ(A_⊥)`.

**Scope.** With governors or `D ≠ 0`, only `span{R_x}` remains.

### A3. Minimal incompatibility hypergraph (TX4 Def. 2, Prop. 1–2; claim V02)

    U(θ, ρ) = {S : α_⊥(S; θ, ρ) ≥ 0},
    H(θ, ρ) = inclusion-minimal elements of U,
    κ = min_{e ∈ H} |e|.

- `H` is an antichain.
- `U` need not be upward closed.
- With a stable base, a target `T` is safe in every implementation order iff it
  contains no hyperedge.

### A4. Boundary localization (TX4 Thm 2)

`H` can change only when some `A_⊥(S)` has an eigenvalue on the imaginary axis.

### A5. Descriptor-affine local actions and exact network closure (TX4 Lemma 1, Thm 3, Cor. 1; claim V08)

The port quantities are

    T_S(s) = g_z + g_x (sI − f_x)^{-1} f_z,
    K(s)   = Eᵀ T_0(s)^{-1} E,
    M      = D K,
    Q      = (I + D K_d)^{-1} D K_o,

and they satisfy

    det P_S = h_S det T_0 det(I + M_SS),
    det(I + M_SS) = Π_i det(I + M_ii) · det(I + Q_SS).

Here `μ_T` is the principal-minor Möbius coefficient, and boundaries are zeros of
`det(I + Q_SS)`.

### A6. Boundary sensitivity (TX4 Prop. 3, eq:grad; claim V10)

For a simple critical eigenvalue,

    ∂s*/∂θ = − p^H ∂_θ M_SS q / p^H ∂_s M_SS q.

**Equivalent forms.** Two forms give the same number:
- the zero form: `det T_S(s*) = 0`, `T_S q = 0`, `p^H T_S = 0`, with
  `∂s*/∂a = − p^H ∂_a T_S q / p^H ∂_s T_S q`;
- the eigenvalue form: `dλ = v^H dA u / v^H u`.

This is the Schur-complement identity `det T_S(s) = det g_z det(sI − A)/det(sI − f_x)`,
valid away from the device poles.

### A7. Symmetry-deflated zero-frequency port (TX4 Thm 4; claim T08)

Complexity: NP-completeness of minimum-cardinality destabilization (claim V26).

### A8. Connected cumulants (claim V22)

These are mathematically valid. They stay **secondary** in CDW, and no bridge or
criticality notion is built on them.

### A9. Classical results used as tools (credited, not claimed)

- implicit function theorem;
- simple-eigenvalue perturbation;
- Gordan's alternative (a common descent direction exists iff `0 ∉ conv{g_j}`);
- Rouché / Gohberg–Sigal argument principle for analytic matrix functions;
- Shapley value;
- submodularity definitions.

---

## B. Inherited evidence (IEEE-39 frozen model; empirical and benchmark-specific)

Source: `results/20260911_FINAL_VALIDATION_MATRIX.csv`. The statuses are
unchanged here.

| TX4 id | evidence | status | what it does NOT prove |
|---|---|---|---|
| V03 | At P4 all 15 proper subsets of {30,33,35,37} are stable and H4 is unstable (α = +0.1270) | VALIDATED, model-specific witness | any general weak-bus set |
| V06, V07 | at fixed network and dispatch, the Q/V policy alone changes H and κ | VALIDATED | that static metrics are "wrong" in general |
| V11 | the port derivative ranks 12 frozen lines' finite ×1.5 effect (Spearman 0.937, one condition); finite sign 9/12 | VALIDATED, local | ranking at other policies or envelopes (**CDW tests this**) |
| V12 | finite-step magnitudes are nonlinear (Spearman 0.2) | NEGATIVE | — |
| V13 | remediation families at P4; no single line ≤ 2× restores P4 | VALIDATED, model-specific | topology impotence in general |
| V14, V29 | governors stabilize P4; κ = 4 is not robust to primary control | BENCHMARK-SPECIFIC / REFUTED | — |
| V15, V16 | non-composability persists (R1 0.62/0.88/1.00/0.84); the witness does not (R2 0.34/0.69/1.00/0.62) | PASS (3/4) / NOT ROBUST | robust contextuality (**CDW tests this**) |
| V28 | the GFL branch ranking does not transfer to the ANDES library GFL (7/12) | REFUTED (transfer) | — |
| V27, N01, N05, N06 | static quantities, simple cycles, effective rank and dominant paths do not explain failures | REFUTED | — |
| dynamic-forest report | frozen-equilibrium line derivative ρ = 0.615 vs re-equilibrated ρ = 1.000 on 12 lines | VALIDATED (one condition) | the frozen/total gap elsewhere (**CDW tests this**) |

**Motivating observation (never run).** The ParaEMT EMT06 preregistration
(`docs/20260911_PAREMT_EMT_PREREG_V1.md:363–368`) wrote the nominal P4 pattern
"each replacement stabilizing alone, destabilizing last" as a hypothesis. It is
a nominal phasor pattern at one policy point and is **not** evidence of
recurrence. No TX4 document tests submodularity.

---

## C. New definitions (CDW)

Throughout:
- the context `S ⊆ V` is the current portfolio;
- the intervention `i ∉ S` is the replacement of SG `i` by the frozen GFL;
- `α(S) = α_⊥(S; θ, ρ, π)`.

### C1. Contextual marginal effect

    Δ_i α(S) = α(S ∪ {i}) − α(S).

- Negative means that replacing `i` in context `S` stabilizes (lowers the
  spectral abscissa); positive means that it destabilizes.
- **Sign classes**, with material threshold `τ_mat` (preregistered):
  - STAB: `Δ ≤ −τ_mat`;
  - DESTAB: `Δ ≥ τ_mat`;
  - NEUTRAL: otherwise.

### C2. Contextual sign reversal

Intervention `i` reverses at `(θ, ρ, π)` if there are contexts `S1`, `S2` with
`Δ_i α(S1) ≤ −τ_mat` and `Δ_i α(S2) ≥ τ_mat`.

- **Stable-context reversal:** `S1` and `S2` are both transversely stable.
- **Global vs mode-tracked reversal.** `α` is a maximum over modes, so a sign
  change of `Δ_i α` can come from a switch of the rightmost mode.
  - Let `λ_c(S)` be the rightmost transverse eigenvalue of `S`.
  - The **tracked marginal** is `Δ_i^tr α(S) = Re λ_m(S ∪ {i}) − Re λ_c(S)`,
    where `λ_m` is the eigenvalue of `S ∪ {i}` matched to `λ_c(S)`.
  - Matching uses the bus-voltage mode shape `φ = −g_z^{-1} g_x u`, which is a
    common coordinate for every portfolio: MAC(φ_S, φ_{S∪i}) ≥ MAC_min and
    `|Δf| ≤ Δf_max`.
  - **Mode-tracked reversal:** a reversal of `Δ_i^tr` in which every
    participating pair is matched and the matched mode is the rightmost mode of
    `S ∪ {i}` (a same-mode transition).

### C3. Weakness objects

These are distinct mathematical objects, and their coincidence is **tested, not
assumed**.

- **W_static** (network and operating point only):
  - Thevenin |Z|, SCR, gSCR;
  - `dV/dQ` at the bus (from the PF Jacobian);
  - electrical distance, effective resistance;
  - Fiedler entry and edge score;
  - edge betweenness;
  - branch |P| and |S| flow, |z_e|.
- **W_dynamic.**
  - For a node or control parameter `a`: `W_a = d Re s*/da`, keeping the sign
    and in both variants of C4.
  - For a branch strength `γ_e`: `W_e = d Re s*/dγ_e`.
  - Normalized per admissible range: `W̄_a = W_a · (a_max − a_min)`.
- **W_contextual:** the distribution `{Δ_i α(S)}_S` and its sign classes,
  variance, sign entropy and reversal indicator. It is a property of the
  pair (intervention, context), **not** of a bus.
- **W_closure** (derived from the exact closure, not from centrality).
  - At a closure zero, `ds*/da = −(∂_a F_Q)/(∂_s F_Q)` with
    `F_Q = det(I + Q_SS)`.
  - With `Q = (I + D K_d)^{-1} D K_o`, the product rule splits `∂_a Q` exactly
    into three parts:
    - `∂D` (device and operating point);
    - `∂K_d` (driving-point coupling);
    - `∂K_o` (transfer coupling among the replaced buses).
  - The **transfer-coupling share** of a branch is the `∂K_o` part of
    `d Re s*/dγ_e`.
  - Physical ablations allowed: only admissible network actions (branch
    reinforcement or outage with PF feasibility). Setting blocks of `K_o` to
    zero is not a physical network and is not used as an ablation.
- **W_control:** the leverage of an admissible parameter, i.e. the finite change
  `Δ Re s*` achievable within its preregistered admissible range. It is
  validated by finite re-solves, and by Newton or continuation where a boundary
  is involved.

### C4. Frozen-operating-point vs re-equilibrated total sensitivity (mandatory distinction)

**Operating point.** Let `w = (x, z, r)` collect the device states, bus voltages
and device references (setpoints) that are solved at equilibrium. Equilibrium
is `R(w, a) = 0`. Two physically distinct semantics are preregistered.

- **SPR (schedule-preserving re-initialization; the TX4 semantics).**
  - Every generator bus keeps its scheduled `P` and `|V|`, and the slack keeps
    `|V|` and its angle.
  - Every device reference is re-derived so that the new operating point is an
    equilibrium. This is exactly what `solve_case` does.
  - Unknowns are `(x, z, r_all)`.
  - Equations are `f = 0`, `g = 0`, plus the schedule equations:
    - SG: `Re S = P_sched` and `|V| = V_set`;
    - GFL: `Re S = P_sched`, `|V| = V_set` and `v_ref = |V|`;
    - slack: `|V| = V_s` and `∠V = θ_s`.
- **RP (reference-preserving).**
  - All device references stay at their base values (`p_m`, `v_ref`, `p_ref`,
    `q_ref`, GFL `v_ref`).
  - Only the slack machine's `p_m` is free, to close the power balance, and one
    gauge equation (`∠V_slack = θ_s`) removes the rotation symmetry.
  - This is the operating point reached if nobody resets setpoints after the
    change.

In both semantics:

    dw*/da = − R_w^{-1} R_a                       (IFT),
    (dA/da)_frozen = ∂_a A |_{w = w*(a0)},
    (dA/da)_total  = ∂_a A + D_w A [dw*/da],
    dλ/da = v^H (dA/da) u / v^H u.

The T-form is `ds*/da = − p^H [T_a + D_w T(dw*/da)] q / p^H T_s q`.

**Structural remarks.** These follow from the model equations. They are exact
and are checked numerically in E3; they are not IEEE-39 findings.

- **(R1)** Under SPR, a controller parameter that enters only dynamic terms
  does not move `(x*, z*)`: PLL gains, PI gains, washout and lead-lag
  constants, the Q/V gain `g_i` (since `e = 0` at the initialized point), and
  `K_A`. It moves only references that do not enter the Jacobian (`v_ref`,
  `p_ref`). Therefore **frozen = total** for these coordinates under SPR.
- **(R2)** A setpoint (`V_set`, `v_ref`, `q_ref`) enters `f` or `g` affinely
  through a constant. It does not enter the Jacobian, so the **frozen partial is
  exactly 0** and all of its effect is re-equilibration.
- **(R3)** Under RP, `K_A` moves the equilibrium, because
  `e_fd = K_A (v_ref − |V|)` at fixed `v_ref`. Network, load and setpoint
  coordinates move the equilibrium under both semantics.

The frozen/total question is therefore substantive exactly for:
- network and topology coordinates;
- operating-point coordinates;
- proportional steady-state controls under RP.

### C5. Plan-level margin

    Φ_T(a) = max_{S ⊆ T} α_⊥(S; a).

- `Φ_T < 0` iff the plan is any-order safe (A3).
- The active set is `J_ε = {S : α(S) ≥ Φ_T − ε}`.
- **Local simultaneous improvement** exists iff some `d` satisfies `g_jᵀ d < 0`
  for all active `j` (in the tangent space of the admissible set).
- **Local conflict witness:** `λ ≥ 0`, `Σ λ_j = 1`, `Σ λ_j g_j = 0` (Gordan).
  It is local only and never a global impossibility claim.

### C6. Local stability-equivalent exchange rate

At a point where `α` of a fixed target is differentiable,

    r_{e←i} = dγ_e/dθ_i |_{α const} = −(dα/dθ_i)/(dα/dγ_e),

using **total** derivatives. It is a local equivalence, not an economic rate.

### C7. Submodularity on the portfolio lattice

- **Submodular:** `Δ_i α(S) ≥ Δ_i α(T)` for all `S ⊆ T`, `i ∉ T`.
- **Supermodular:** the reverse inequality.
- A **violation** is a triple `(S, T, i)` that breaks the inequality by at least
  `τ_res`.
- A **minimal counterexample** has `|T ∖ S| = 1`. It suffices to check these:
  a set function is submodular iff the one-step condition
  `Δ_i α(S) ≥ Δ_i α(S ∪ {j})` holds for all `i ≠ j ∉ S`.

### C8. Graph-mode operators (only where physically defined)

`B_bus` is decomposed exactly as

    B_bus = L_B + Δ_tap + B_charging + B_shunt,

where:
- `L_B` is the Laplacian of branch coupling weights `b_e/t_e`;
- `Δ_tap` is the exact diagonal residual from off-nominal taps;
- `B_charging` and `B_shunt` are diagonal.

The graph Fourier basis is `U` from `L_B = U Λ Uᵀ`. The mixing object is

    T̂(s) = (Uᵀ ⊗ I₂) T(s) (U ⊗ I₂),
    R_mix = T̂ − blockdiag(T̂_kk),
    μ_mix = ‖R_mix‖_F / ‖T̂‖_F.

### C9. Decision-preserving reduction

- A reduced model is admissible only if it:
  - has materially fewer states (≤ 70 % of the full dimension);
  - does not evaluate the eliminated full complement online.
- It is judged on verdicts, H, κ and intervention rankings on holdout data.
- **Certificate.** With `R_S = I + D_S K_SS` and `R̃_S` its reduced counterpart,
  the condition `sup_Γ ‖R̃_S^{-1}(R_S − R̃_S)‖ < 1` on a contour `Γ` enclosing no
  poles of either function gives the same number of zeros inside `Γ`
  (Gohberg–Sigal).
  - A **sampled** check of this condition is labelled "SAMPLED, NON-RIGOROUS".
    Otherwise the answer is ABSTAIN.

---

## D. New hypotheses (to be tested; see `docs/CDW_PREREG_V1.md`)

| id | hypothesis |
|---|---|
| H1 | Contextual sign reversal: some interventions change sign across legitimate contexts, recurrently across policies, and it survives mode tracking |
| H1b | The reversal is strong enough to make every fixed node-only ranking materially wrong for next-replacement decisions |
| H1c | Contextuality persists under parameter envelopes even where the minimal failing set changes (robust phenomenon ≠ robust identity) |
| H2 | Total re-equilibrated line and node sensitivities differ materially from frozen partials, and they predict held-out finite effects better than static baselines |
| H3 | Robust dynamic weak corridors exist: admissible corridor coordinates whose ranking persists across holdout policies and envelopes |
| H4 | Topology actions alone can remove or create incompatibilities, and static topology scores do not predict which |
| H5 | Plan-level multi-constraint design succeeds where single-boundary tuning leaves an unsafe intermediate subset |
| H6 | A low-dimensional graph-mode or dynamic basis preserves decisions (spectral closure); the policy changes modal mixing at fixed L |
| H7 | A genuinely reduced model preserves decisions with speedup and a certificate or abstention rule |
| H8 | Dominant observed modal energy is often not the stability-limiting mode (phasor-domain) |
| H9 | Findings transfer to an alternative converter model |
| HS | α_⊥ is neither submodular nor supermodular on the tested class; a tracked-mode margin is no more structured |

**Central meta-hypothesis to falsify:** static weakness ≠ dynamic weakness ≠
contextual weakness ≠ intervention leverage.

---

## E. Not yet tested (all CDW results)

Every hypothesis in D has status HYPOTHESIS or NOT_YET_TESTED until the
preregistered experiments run. Explicitly **not** claimed at any stage:
- universal weak buses or lines;
- probabilities from envelope coverage;
- static topology as sufficient;
- simple-cycle or connected-cumulant causality;
- transfer across converter models without a holdout;
- that a better `α_⊥` implies a larger nonlinear region of attraction (it does
  not, mathematically; see E22 placeholder);
- EMT validity (the phasor-domain model only).
