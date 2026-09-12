# Author decision memo: three open manuscript choices (Phase 11)

Date: 2026-09-11.

- Target manuscript:
  `reports/papers/tx4_policy_dependent_incompatibility/main.tex` (TX4, a
  single file).
- **The manuscript has not been edited.** The diffs below are proposals. They
  are to be applied only after the author decides, and only now that the final
  validation ledger (`docs/20260911_FINAL_VALIDATION_LEDGER.md`) is complete.

---

## A. The pandapower fix — recommendation: **ACCEPT FIX**

### What happened

The supplied cross-tool pack converts the canonical JSON to a PYPOWER branch
table (`validation_runs/cross_tool_final/pack/canonical_ppc.py:71`, copying
`tap`). pandapower then builds the network with `from_ppc`
(`pack/run_pandapower_parity.py:23`).

For a branch whose from-bus is the LV bus, `from_ppc` swaps `hv_bus`/`lv_bus`
and sets `tap_side = "hv"`. The tap therefore ends up on the other terminal
(`from_ppc.py:229–238, 253, 285–292` in pandapower 3.4.0).

Three branches are affected; all three have b = g = 0 and φ = 0:

| idx | canonical (from→to, r, x, t) | translated (from→to, r′, x′, t′) |
|---|---|---|
| 35 | 31→6, 0, 0.025, 0.9 | 6→31, 0, 0.02025, 1/0.9 = 1.1111 |
| 37 | 12→11, 0.0016, 0.0435, 1.006 | 11→12, 0.0016193, 0.0440236, 0.99404 |
| 38 | 12→13, 0.0016, 0.0435, 1.006 | 13→12, 0.0016193, 0.0440236, 0.99404 |

### Transformer equations

Let `y = 1/(r + jx)`, and let the off-nominal tap `t` sit on the from bus f
(canonical, MATPOWER, ANDES `Line.build_ybus` and the project model):

```
canonical (tap t at f):          [ I_f ]   [  y/t²   −y/t ] [ V_f ]
                                 [ I_t ] = [ −y/t     y   ] [ V_t ]

as supplied (tap moved to t):    [ I_f ]   [   y     −y/t ] [ V_f ]
                                 [ I_t ] = [ −y/t    y/t² ] [ V_t ]      ← a different two-port

translated (HV bus as from,      y' = y/t²,  t' = 1/t, from = t-bus:
tap t' at the HV bus):           Y'_tt = y'/t'² = y,  Y'_tf = −y'/t' = −y/t,  Y'_ff = y' = y/t²
                                 ≡ canonical, exactly
```

- **Numerical signature of the as-supplied error.** For branch 35,
  `y = −j40`. The as-supplied pack puts `y/t² = −j49.383` at bus 6 instead of
  `−j40`: a 9.383 pu Ybus error at (6,6) and (31,31), which is exactly what
  was measured. The measured error at (12,12) is 0.55 pu.
- **After translation.** The Ybus agrees with the project model to 1.14e-13,
  and the power flow to 3e-14 in Vm/Va. Branch flows agree to 2e-10 MVA
  (`docs/FINAL_CROSS_TOOL_VALIDATION.md` lines 57–92;
  `results/PANDAPOWER_STATIC_PARITY.csv`).
- **Script checks.** `validation_runs/cross_tool_final/run_pack_pandapower_tapside.py`
  asserts two-port equality before use; the assert gate is 1e-12, and the
  measured value is 1.4e-14.

### Classification

**An equation-preserving convention translation, not a model adjustment.**
- The translated rows are an exact algebraic re-expression of the same ideal
  transformer and series impedance, and no parameter was fitted.
- The *as-supplied* conversion is the one that changed the physics: it moved
  the tap to the other winding without re-referring the impedance.
- The three other compatibility fixes (import path, `p/v/sn` key aliases, and
  the pandapower 3.4.0 RATE_A = 0 indexing bug) are non-electrical.
- The Ybus check confirms that none of them changes the network.

**Justification for ACCEPT.** The Ybus identity to 1e-13 makes the fix
verifiable by anyone. The alternative, rejecting it, discards a correct
independent static reproduction because of a converter-script convention.

**Required wording if accepted.** The manuscript must say that the
pandapower representation required a tap-convention translation. The
translation is included in diff B below.

**If rejected.** Use only the Python–ANDES static agreement: delete the
bracketed pandapower clause in diff B.

---

## B. Limitations sentence — the current "within 4.1 %" is stale

### Current text (`main.tex:959–962`)

```latex
\item \emph{Independent implementation.} ANDES reproduces the power flow, the
base inter-area mode (within 4.1\%), and the governor damping direction. It
does not implement the converter model, so no portfolio result is claimed as
independently reproduced.
```

### Why it is stale

Two campaigns have replaced it:

1. The final equation-equivalent cross-tool campaign (0a4e18c2):
   - identical Ybus to 1e-13;
   - bus voltages within 1.7e-7 (ANDES's internal +1e-8 r/x regularization);
   - R3 electromechanical eigenvalues within 1.5e-6 s⁻¹ on 12 subsets;
   - branch sensitivities 12/12 signs, Spearman 1.00;
   - port derivative 11/12, Spearman 1.00.
2. The Phase 8 custom-GFL reproduction (7773ea8c):
   - the project's own two-axis machine and GFL equations were transcribed
     into ANDES, identical to 1.5e-14 at device level;
   - all 32 verdicts reproduced (16 subsets of H4 at P4 and at g = 0.25);
   - H = {H4} and κ = 4 at P4, H = ∅ at g = 0.25;
   - max |Δα⊥| 1.3e-6 s⁻¹.

The sentence "no portfolio result is claimed as independently reproduced" is
now too weak. A sentence claiming "ANDES validates the converter model" would be
too strong.

### Proposed diff (accept-fix version)

```diff
-\item \emph{Independent implementation.} ANDES reproduces the power flow, the
-base inter-area mode (within 4.1\%), and the governor damping direction. It
-does not implement the converter model, so no portfolio result is claimed as
-independently reproduced.
+\item \emph{Independent implementation.} The network and operating point were
+reproduced in ANDES~2.0.0 and, after a tap-convention translation of the three
+transformers whose off-nominal tap is on the low-voltage from-bus, in
+pandapower~3.4.0 (bus-admittance matrices identical to $10^{-13}$~p.u.; bus
+voltages within $2\times10^{-7}$~p.u., the residual being ANDES's internal
+$10^{-8}$~p.u. branch-impedance regularization). On an equation-equivalent
+synchronous configuration ANDES reproduces the electromechanical eigenvalues
+within $2\times10^{-6}$~s$^{-1}$ and the branch sensitivities of the critical
+mode on a preregistered 12-branch holdout (12/12 signs, Spearman 1.00). With
+the two-axis machine and the grid-following converter of this paper
+transcribed into ANDES as custom models, ANDES reproduces all 32 transverse
+verdicts of the 16 subsets of $\{30,33,35,37\}$ at P4 and at $g=0.25$
+($|\Delta\alpha_\perp|\le1.3\times10^{-6}$~s$^{-1}$), including
+$\Hy=\{\{30,33,35,37\}\}$ and $\kappa=4$ at P4. This independently reproduces
+the computation, not the adequacy of the converter model: no conclusion is
+claimed for other converter models, and with a library converter model the
+branch ranking changes (7/12 signs). ANDES also confirms the governor damping
+direction.
```

**Reject-fix version.** Replace the first sentence with: "The network and
operating point were reproduced in ANDES~2.0.0 (bus-admittance matrices
identical to $10^{-13}$~p.u.; …)". That is, drop the pandapower clause and
keep the rest.

**Wording that must not be used:**
- "ANDES independently validates the GFL model."
- "validated against a second converter model."
- Any claim for P_inf. P_inf was optional in Phase 8 and was not run.

---

## C. C2 clarification — the full-order term is not an irreducible four-device interaction

### Facts

- **Boolean statement at s\*.** The statement "no interaction truncation
  below the cardinality of the changing minimal coalition reaches the
  characteristic zero" concerns the Boolean (Möbius) truncations **evaluated at
  the exact boundary zero s\***. It remains exactly true (main.tex:643–655;
  frozen FC18).
- **Pair identity.** Because Q has zero 2×2 diagonal blocks, the full-order
  Boolean coefficient of a four-unit set is

  `μ_{1234} = χ_{1234} + χ_{12}χ_{34} + χ_{13}χ_{24} + χ_{14}χ_{23}`,

  i.e. the connected four-way cumulant plus the products of pair cumulants
  that jointly cover all four buses. This was checked numerically to 3e-14 on
  a random Q with zero diagonal blocks.
- **The flagship is composite.** The connected-cumulant audit
  (docs/CONNECTED_CUMULANTS_EXTENSION.md, FLAG row and lines 76–90; claim C07)
  classifies the flagship P4-line
  boundary as COMPOSITE:
  - ν = 0.069 < 0.10;
  - deleting χ_{30,33,35,37} moves the zero by 0.024 s⁻¹;
  - the dominant connected cluster is {33, 35}.

  The label is threshold-sensitive: the boundary would be CONNECTED at 0.05.
- **Lower orders at P4.** At P4 itself (PCV03), the zero of the order-≤3
  Boolean truncation of det(I + Q) lies at Re s = +0.1267, against the exact
  +0.1270. The order-≤2 zero is at −0.0148, and every spectral-abscissa
  truncation predicts stability:
  - B3 −0.203;
  - B4 −0.213;
  - B5 −0.210;
  - B5b −0.184.

  So "only the full order closes the determinant at s\*" must not be read as
  "no lower-order computation can flag H4".

### Proposed diffs

**C2 contribution (`main.tex:125–127`):**

```diff
-(Theorem~\ref{thm:factorization}, Corollary~\ref{cor:closure}). At all ten
-frozen IEEE-39 boundaries, only the full-order interaction term of the changing
-coalition closes the determinant.
+(Theorem~\ref{thm:factorization}, Corollary~\ref{cor:closure}). At all ten
+frozen IEEE-39 boundaries, only the full-order interaction term of the changing
+coalition closes the determinant at the boundary zero; this full-order term is
+a closure order, not an irreducible $|S|$-device interaction.
```

**Section IV-B, insert after `main.tex:655`** ("…benchmark."):

```diff
 $(\theta,s)$ draws (suite P3), the full-order coefficient never vanishes
 ($7.8\times10^{-4}$ to 469), so no lower-degree representation is exact on this
 benchmark.
+
+The full-order term must not be read as an irreducible four-device
+interaction: the minimum destabilizing cardinality $\kappa$ and the order of the
+irreducible (connected) interaction are different quantities. Because $Q$ has
+zero diagonal blocks, the four-unit coefficient is
+$\mu_{1234}=\chi_{1234}+\chi_{12}\chi_{34}+\chi_{13}\chi_{24}+\chi_{14}\chi_{23}$,
+where $\chi$ are the connected (partition-lattice) cumulants of
+$\det(I+Q_{SS})$. At the four-unit boundary the connected four-device term
+carries a normalized weight of only $0.069$, and deleting it moves the
+characteristic zero by $0.024$~s$^{-1}$; the dominant connected cluster is
+$\{33,35\}$. The four-unit boundary is therefore composite: it is closed by
+products of lower-order connected interactions that jointly touch all four
+buses. For the same reason, at P4 the zero of the order-$\le3$ truncation lies
+at $+0.1267$~s$^{-1}$, close to the exact $+0.1270$~s$^{-1}$, whereas
+additive, pairwise and third-order truncations of the spectral abscissa all
+predict stability ($-0.18$ to $-0.21$~s$^{-1}$).
```

**Figure caption (`main.tex:690–691`):**

```diff
-interaction expansion \eqref{eq:mu}. Only the full-order term (order $|S|$)
-reaches the characteristic zero. For pairs, order $\le2$ is the full order.}
+interaction expansion \eqref{eq:mu}. Only the full-order term (order $|S|$)
+reaches the characteristic zero; at the four-unit boundary this term is
+dominated by products of pair interactions, not by an irreducible four-device
+term. For pairs, order $\le2$ is the full order.}
```

The abstract (lines 67–69) and the conclusion (lines 1021–1023) are correct as
Boolean statements and need no change. Optionally, add "(a closure order, not
an irreducible interaction order)" after "cardinality of the changing minimal
coalition" in the conclusion.

---

## Summary of recommendations

| item | recommendation | blocking? |
|---|---|---|
| A. pandapower | **ACCEPT FIX** (equation-preserving tap-convention translation; Ybus identity 1e-13); state it in the text | no |
| B. Limitations | **REPLACE** with diff B (accept-fix version); it now includes the Phase 8 custom-GFL reproduction, scoped to reproducing the computation | no, but the current text is factually stale |
| C. C2 | **ADD** the clarification (diffs C); keep the Boolean statement | no |
