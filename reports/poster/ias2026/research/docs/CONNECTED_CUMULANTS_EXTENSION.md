# Finite-amplitude connected cumulants: report (Phases 1–7)

**Governing documents.**

- The theory is `theory/FINITE_AMPLITUDE_CONNECTED_CUMULANTS.md`.
- The preregistration is `configs/ias2026/connected_cumulants_prereg_v1.yaml`,
  committed as de09db3b **before** any IEEE-39 cumulant was computed.
- The canonical ledger is `docs/20260911_CANONICAL_THEORY_LEDGER.md`, with its
  matrix in `results/20260911_CLAIM_MATRIX.csv`.

**Frozen-data rules.**

- Frozen data only: the ten FC18 boundaries (θ and s* read verbatim), the FC04
  Pg-matched set, the 418 frozen points of the F7B line, and the 12 frozen
  holdout lines.
- No new search and no retuning.

## Phase 1: mathematics (proved)

| item | result |
|---|---|
| C1 | Moment–cumulant inversion `F(S) = Σ_π Π chi_B` and uniqueness, via the exponential formula in `C[x]/(x_i²)` (classical) |
| C2 | Factorization across a partition ⇔ **all** mixed cumulants vanish. The top cumulant vanishing alone does not imply factorization (path counterexample). |
| new, elementary | A minimal unstable coalition is **cumulant-connected** at its characteristic zero, a pointwise algebraic statement. At a boundary, `chi_H(s*) = −C_H(s*)`, so `|chi_H|` there is forced and uninformative by itself. |
| C3, scalar | **TRUE** for `|S| ≥ 2`: `chi_S = (−1)^{|S|−1} Σ spanning directed cycles`. For `|S| = 1`, `chi_{i} = 1 + q_ii`, not the cycle value, so the singleton case is excluded. It coincides with the multilinear coefficient of `log det(I + ΘQ)`. |
| C3, 2×2 ports | The scalar formula does **not** transfer. The exact statement is the block-connected permutation expansion, e.g. `chi_ab = −tr(Q_ab Q_ba) + det Q_ab · det Q_ba`. The naive holonomy sum is only the infinitesimal cumulant. |
| C4 | `d_alg`, `d_conn` and κ are logically independent. Five constructions realize every non-implication; for example κ = 4 with `d_conn = 2` (path) and κ = 4 with `d_alg = 2` (star). |

## Phase 2: toy falsification (`results/CC/CC02`, 74 unit tests)

All gated residuals are at machine precision:

| quantity | residual |
|---|---|
| C1 inversion | ≤ 1.9e-13 |
| C2 mixed cumulants | ≤ 5.3e-15 |
| scalar cycle formula (n ≤ 6) | ≤ 6.2e-14 |
| block expansion | ≤ 5.3e-15 |
| two-block closed form | ≤ 7.8e-16 |
| gauge drift | ≤ 7.7e-15 |

The naive holonomy formula for 2×2 ports is off by ≥ 4.5 % (relative).

- **T1.** Two independent pairs: `|μ_1234| ≥ 5e-3` while `|chi_1234| ≤ 6e-17`.
- **T2.** Ring: `chi_1234 = −q12 q23 q34 q41` to 3e-16.
- **T3.** Same κ = 4, two structures: the symmetric path is **composite**
  (ν = 0, `chi_V ≡ 0`) and the directed ring is **connected** (ν = 1/2).
- **Connectivity theorem.** Random toys give 184 minimal coalitions, with 0
  violations; 5 of them are composite.

## Phase 3: the ten frozen boundaries (`results/CC/CC03`)

- The reconstruction `F(H) = Σ_π Π chi_B` holds to ≤ 1e-16, and
  `|F(H)(s*)| ≤ 4e-13`.
- Theorem 2.3 holds at 10/10 boundaries.
- The C2 and PortActionSpace machineries agree on `F` to 4e-9.

| boundary | H | ν_H | label at 0.10 (primary) | at 0.05 / 0.20 | dominant δ cluster | deletion displacement (s⁻¹) | materially necessary |
|---|---|---|---|---|---|---|---|
| E1 κ 4→3 | 30+33+35 | 0.147 | CONNECTED | C / composite | 30+33 | 0.072 | yes |
| E2 | 30+33+37 | 0.183 | CONNECTED | C / composite | 30+33 | 0.064 | yes |
| E3 | 33+35+37 | 0.182 | CONNECTED | C / composite | 33+35 | 0.038 | no |
| E4 κ 3→2 | 30+33 | 0.5 | forced (pair) | — | — | no zero | — |
| E5 | 33+35+37 | 0.147 | CONNECTED | C / composite | 33+35 | 0.039 | no |
| E6 κ 2→3 | 30+33 | 0.5 | forced (pair) | — | — | no zero | — |
| E7 | 30+33+35 | 0.183 | CONNECTED | C / composite | 33+35 | 0.056 | yes |
| E8 κ 3→4 | 30+33+37 | 0.141 | CONNECTED | C / composite | 30+33 | 0.029 | no |
| E9 → ∅ | 30+33+35+37 | **0.057** | **COMPOSITE** | C / composite | 30+33+35 | 0.019 | no |
| FLAG (P4 line, 0.706 Hz) | 30+33+35+37 | **0.069** | **COMPOSITE** | C / composite | 33+35 | 0.024 | **no** |

In the "at 0.05 / 0.20" column, C means CONNECTED.

**Central question.** Is `chi_{30,33,35,37}(jω*)` materially necessary at the
flagship boundary? **No, by the preregistered rule.**

- `ν = 0.069 < 0.10`.
- Deleting the four-way cumulant moves the characteristic zero by 0.024 s⁻¹, to
  Re −0.0055, below the preregistered 0.05 s⁻¹.
- The four-way connected term is non-zero, and the support is connected as the
  theorem requires. But the zero is assembled mainly from lower-order connected
  clusters: `delta` is 0.037 for 33+35, 0.032 for 30+33+35, and 0.022 for the
  full set.

The label is threshold-sensitive: the flagship would be CONNECTED at 0.05.

**Relation to the manuscript's C2 statement.** The statement "no truncation below
order |S| reaches the zero" is a *Boolean* statement. It remains exactly true:
the order-3 truncation leaves `|μ_4| = 0.16`. The cumulant layer shows that this
full-order Boolean term is not an irreducible four-way interaction, and that its
eigenvalue weight is 0.024 s⁻¹.

**Remainder.** At the eight boundaries with `|H| ≥ 3`, the naive single-cycle
(holonomy) part misses 39–65 % of the connected interaction.

## Phase 4: does chi add information? (`results/CC/CC04`)

The data are 10 failing and 25 Pg-matched stable size-4 portfolios at the E12
census policy. Every quantity is evaluated at `j Im λ_P`. The orientations are
preregistered.

| diagnostic | AUC | permutation p |
|---|---|---|
| **`|chi_P|`** | **0.952** | 1e-4 |
| closure distance (existing return-difference test) | 0.944 | 1e-4 |
| `|μ_P|` (full-order Boolean) | 0.812 | 1.5e-3 |
| MVA | 0.716 | 0.026 |
| gSCR | 0.676 | 0.057 |
| ν at the mode | 0.600 | 0.19 |
| min SCR | 0.474 | 0.68 |
| Pg | 0.112 | 1.0 |
| largest simple-cycle holonomy (E18) | **0.084** (anti-predictive) | 1.0 |

**Preregistered verdict: PROMOTE.** AUC 0.952 ≥ 0.80, it exceeds every
non-definitional diagnostic by ≥ 0.05 (the best is 0.812), and p ≤ 0.01.

**Qualification.** Over the whole set, `|chi|` only *matches* the existing
closure distance (0.952 against 0.944).

**Post hoc, descriptive, not preregistered.** The failing portfolios' critical
modes lie at 0.42–0.57 Hz, while most controls' lie at 0.7–1.4 Hz. Within the
0.35–0.60 Hz window (10 failing against 6 stable), `|chi|` reaches AUC 0.983
while the closure distance reaches 0.767. The separation is therefore not a pure
frequency artifact, but n is small. The script is
`experiments/connected_cumulants/CC04b_posthoc_frequency_window.py`, and its
output is `results/CC/CC04/CC04_posthoc_frequency_window.txt`.

**Clarification of the old failure.** Single-cycle holonomies are *larger* for
the stable portfolios. The connected cumulant is a signed sum over all
block-connected permutations, and it points the other way.

## Phase 5: policy mechanism track (`results/CC/CC05`)

The data are all 418 frozen points of the F7B line. `D(g)` is the dominant
connected support of the core's critical mode.

- `|D| = κ` at 73 %.
- `D ∈ H` at **34 %**.
- `D` is 33+35 at 234 points, 30+33+35 at 156, and 30+33 at 10. The four-way
  cluster is never dominant.

**H5 is REFUTED.** Witness contraction does not correspond to a change of the
dominant connected support.

## Phase 6: diagnosis and design (`results/CC/CC06`)

The tools are `dF(B)/da = tr(adj(I+Q_BB) dQ_BB/da)` and the partition product
rule, which remain valid at the singular boundary.

**Recomputation checks.**

- The recomputed `ds*/da` matches the frozen FC18 values to 2e-6.
- The line `ds*/dγ` matches the frozen F2c predictions: 12/12 signs,
  Spearman 1.0.

**Policy coordinates (10 boundaries).** The top-1 argmins coincide at 9/10, so
the preregistered verdict is SUPPORTED. The agreement is not mechanistic:

- the two pair boundaries agree by necessity (σ ≡ 1);
- the excitation time-constant scale `t` dominates both rankings nearly
  everywhere;
- **at the flagship they disagree.** The most stabilizing coordinate is the Q/V
  gain `g` (Re ds*/dg = −0.46). Increasing `g` *raises* `|chi_V|` (+0.35).

**12 frozen holdout lines (flagship boundary).**

- The most stabilizing line is L35 (31–6). The most `|chi_V|`-suppressing line
  is L17 (13–14). Spearman is −0.34.
- **H6 is REFUTED.**
- Reinforcing 11 of the 12 lines *increases* `|chi_V|`.
- Only a median 16 % of the boundary motion passes through the connected
  four-way term (`|σ|`).

## Phase 7: independent GFL reproduction

Phase 7 was **not executed.** A correct execution needs a separate ANDES copy and
environment, because the frozen `tx3-andes` runs from the vendored tree, plus a
device-level equation-equivalence test before any system comparison. The ledger
§7 gives the details. The status is NOT YET TESTED, and it is not a blocker.

## Deviations and run notes (disclosed)

1. **First CC03 launch.** BLAS thread variables were set after numpy was
   imported, which made solves about 300× slower. I stopped my own worker
   processes (a foreign job was left untouched) and relaunched with
   `OPENBLAS/OMP/MKL_NUM_THREADS=1` set in the shell. No result came from the
   aborted run.
2. **Second launch.** It stopped on `|H| = 2`, where `C_H ≡ 1` has no zero. A
   guard reports "C_H constant (no zero)". The preregistered material-necessity
   rule applies only to `|H| ≥ 3`, so no criterion was affected.
3. **Lazy Jacobians.** `_cc.Realization` evaluates only the C2 vertex Jacobians
   it uses. These are the same calls as FC18; the change is for performance only.
4. **Post hoc analyses** (the frequency window, and the `|chi|`–closure
   correlation of −0.62) are labelled descriptive and play no part in any
   verdict.

## Answers

1. **Mathematically proved:**
   - classical: the transverse quotient; the hypergraph axioms; the planning
     hierarchy; the Metzler class; the complexity corollary; the resilience
     complex; the principal-minor structure; the zero-frequency port; the
     Sylvester ratio identity; the port-derivative formula;
   - new, elementary: C1, C2 and its converse; cumulant-connectivity of minimal
     coalitions and the boundary identity; the scalar cycle formula
     (`|S| ≥ 2`); the block-connected expansion and the refutation of the naive
     holonomy formula; the independence of `d_alg`, `d_conn` and κ.
2. **Toy / synthetic only:**
   - the existence of both composite and connected κ = 4 coalitions (path
     against ring);
   - the C4 separations;
   - the random-toy connectivity checks;
   - the MC01 S1–S7 generality checks;
   - the zero-frequency port holdout (synthetic plus Kundur; no IEEE-39 real
     crossing).
3. **Established on IEEE-39 (benchmark-specific or conditional):**
   - policy-dependent minimal incompatibility and witness contraction;
   - the flagship κ = 4 at P4 and ∅ at P_inf;
   - the governed replication;
   - the inter-area mechanism;
   - the condenser rating;
   - the order-dependent targets;
   - the nonlinear = spectral result inside scope;
   - subcritical Hopf;
   - the network-closure anatomy;
   - the port-derivative accuracy;
   - the new CC results: composite flagship; connected triples at 0.10;
     `|chi|` AUC 0.95; H5 and H6 refuted.
4. **Independently reproduced:**
   - network and equilibrium (ANDES; pandapower after the tap-convention fix);
   - equation-equivalent synchronous dynamics;
   - equation-equivalent branch sensitivities and the port derivative.
5. **Conditional on the custom GFL:** everything portfolio-level:
   - H, κ, P4, P_inf, the policy maps, the witness contraction, the condenser
     and order results;
   - all CC results on IEEE-39;
   - the GFL line ranking. As a transferable conclusion this is refuted: 7/12
     signs across converter models.
6. **The flagship is composite.**
   - By the preregistered rule it is a COMPOSITE four-replacement interaction:
     ν = 0.069, and `chi_V` is not materially necessary (deletion 0.024 s⁻¹).
   - All four units are needed for instability (κ = 4). The cumulant support is
     connected.
   - The characteristic zero is assembled mostly from the pair and triple
     connected clusters (33+35, 30+33+35).
   - The label is threshold-sensitive: CONNECTED at 0.05.
7. **Actionable value is limited.**
   - *Diagnosis:* modest and benchmark-specific. `|chi|` separates failing from
     matched stable portfolios (preregistered PROMOTE, AUC 0.95). It explains
     why cycle scores fail. But it only matches the existing closure test
     overall.
   - *Mechanism tracking (H5) and line design (H6):* refuted.
   - *Policy agreement:* not mechanistic.
   - Recommendation: at most a short methods paragraph or appendix; not a
     headline and not a design tool.
8. **Before TPWRS submission, nothing new is scientifically necessary.**
   Recommended:
   - one clarifying sentence stating that the full-order Boolean closure term is
     not an irreducible four-way interaction. At the flagship it is composite in
     cumulant terms, with an eigenvalue weight of 0.024 s⁻¹;
   - the author decisions left open by the cross-tool campaign: the pandapower
     fix and the Limitations sentence.

   Optional: the independent GFL reproduction (Phase 7) and the
   machine-parameter sensitivity (old gap 2).
