# PD39 255+1 - Final mechanism and generalization validation

Campaign base: commit 902403cf
Branch: research/pd39-255plus1-mechanism-validation
Preregistration commit: 5f28df82

FINAL CASE: C

## 1. Executive verdict

This campaign did execute the requested new four-gate validation. It does not validate the strongest version of the 255+1 story. The frozen discovery pattern remains exact in the original nine conditions, and the nonlinear high-PLL TDS preserves the fragile/robust ordering. However, the preregistered independent central-FD alpha check fails for 72/81 cases; the complete 24-condition census produces 42 H0 blockers and 47 H0.05 blockers; the physical SG-to-GFL homotopy is unavailable; and the exact K,D,Q closure test is BLOCKED by missing model semantics. Under the preregistered decision rule this is FINAL CASE C: the claims must be narrowed and the IAS poster should not be redesigned around a validated universal mechanism.

## 2. Gate A - numerical truth audit

The audit contains 9 portfolios x 9 scenarios x 3 tolerances = 243 tolerance rows and 81 independent rows. All 243 tolerance rows completed. The maximum alpha difference between 1e-10 and 1e-12 was 7.30647986558e-10 s^-1, with no sign or 0.05-classification changes in that tolerance comparison. Descriptor generalized eigenvalues also agree with the reference reduction; the maximum descriptor difference was 4.45388462222e-10 s^-1.

The preregistered independent central finite-difference Jacobian was not a successful independent alpha validation: 72/81 comparisons exceeded 1e-4 s^-1, with maximum difference 0.138325400017 s^-1. There were 0 true-stability sign changes but 67 engineering-margin classification changes. A one-case step probe showed that the full Jacobian norm can be close while the rightmost eigenvalue is highly sensitive; the V8 nominal eigenvector condition was approximately 5.8e4. This explains the numerical behavior but does not erase the preregistered failure.

Gate A: FAIL on the frozen independent-alpha criterion; tolerance and descriptor subchecks PASS.

## 3. Gate B - high-PLL nonlinear TDS

The corrected implementation applies the common +1% P/Q constant-power-factor pulse at deterministic load bus 39 from 1.0 to 1.1 s and integrates to 20 s. Exactly 12 simulations completed, all with sign consistency.

- V8: mean linear alpha 0.112214916; mean nonlinear estimated rate 0.1119606096
- best_7of8: mean linear alpha -0.1383254208; mean nonlinear estimated rate -0.2850800408
- worst_7of8: mean linear alpha -0.1163587056; mean nonlinear estimated rate -0.1140069162

V8 has positive estimated growth in all four high-PLL cases. Both 7/8 controls decay, with the best 7/8 more strongly damped than the worst 7/8. Gate B: PASS for sign and relative-order consistency. The earlier 50 ms implementation is not used as a result.

## 4. Gate C - complete holdout census

The census contains exactly 6144 rows, 6144 unique portfolio-condition keys, 256 portfolios, and 24 conditions. There are 17 failed equilibrium/spectrum attempts; all are retained and are not silently scored as instability.

Condition-wise counts:

- C01: H0=4, H0.05=4, failed=0
- C02: H0=2, H0.05=2, failed=8
- C03: H0=2, H0.05=2, failed=0
- C04: H0=0, H0.05=0, failed=1
- C05: H0=3, H0.05=5, failed=0
- C06: H0=0, H0.05=0, failed=0
- C07: H0=2, H0.05=2, failed=0
- C08: H0=2, H0.05=3, failed=0
- C09: H0=0, H0.05=0, failed=0
- C10: H0=0, H0.05=0, failed=0
- C11: H0=0, H0.05=0, failed=0
- C12: H0=19, H0.05=19, failed=8
- C13: H0=0, H0.05=0, failed=0
- C14: H0=0, H0.05=0, failed=0
- C15: H0=0, H0.05=0, failed=0
- C16: H0=1, H0.05=2, failed=0
- C17: H0=1, H0.05=1, failed=0
- C18: H0=2, H0.05=2, failed=0
- C19: H0=0, H0.05=0, failed=0
- C20: H0=3, H0.05=4, failed=0
- C21: H0=0, H0.05=0, failed=0
- C22: H0=0, H0.05=0, failed=0
- C23: H0=0, H0.05=0, failed=0
- C24: H0=1, H0.05=1, failed=0

Totals are H0=42 and H0.05=47. Nonzero blockers occur in 12/24 conditions for both definitions. V8 appears in H0 only at C24 and appears in H0.05 at none; therefore the exact V8 H0.05 identity does not generalize. Gate C: FAIL_GENERALIZATION, despite the census execution itself passing its row-count integrity checks.

## 5. Gate D - penetration, composition, and mode mechanism

The eight 7/8 predecessors have the following frozen discovery robust margins:

- missing 34: portfolio 30;32;33;35;36;37;38; converted MW 4112.0; m9 0.0858347015467; nominal alpha -0.0967153322641
- missing 30: portfolio 32;33;34;35;36;37;38; converted MW 4370.0; m9 0.126228983466; nominal alpha -0.130545971751
- missing 35: portfolio 30;32;33;34;36;37;38; converted MW 3970.0; m9 0.136356176224; nominal alpha -0.138251947158
- missing 38: portfolio 30;32;33;34;35;36;37; converted MW 3790.0; m9 0.13801971639; nominal alpha -0.138019716734
- missing 32: portfolio 30;33;34;35;36;37;38; converted MW 3970.0; m9 0.138022315453; nominal alpha -0.138022318748
- missing 36: portfolio 30;32;33;34;35;37;38; converted MW 4060.0; m9 0.138057659073; nominal alpha -0.13805765951
- missing 33: portfolio 30;32;34;35;36;37;38; converted MW 3988.0; m9 0.138295566944; nominal alpha -0.138295569249
- missing 37: portfolio 30;32;33;34;35;36;38; converted MW 4080.0; m9 0.138325419498; nominal alpha -0.138325422033

The 6/8 and 7/8 matched comparisons are aggregate-matched: converted MW and IBR MVA differences are zero within each pair, while alpha differences remain nonzero and scenario-dependent. The 90-row table distinguishes frozen-discovery alpha for 6/8 pairs from AD-reference alpha for the newly audited 7/8 pairs; no finite-difference alpha is used for this interpretation.

All 81 newly tracked critical modes are classified as electromechanical/control by the installed modal participation rule. V8 is a low-frequency oscillatory critical mode at roughly 0.32-0.37 Hz in the four high-PLL cases, with positive alpha roughly 0.110-0.114 s^-1. The 7/8 critical modes are generally distinct in the common bus-voltage observable: MAC values range from 0.00219161 to 0.907793 where available, and several 7/8 critical modes are real/near-zero-frequency or have different oscillatory frequencies. The evidence supports a physical modal distinction, but does not establish a continuous SG-to-GFL homotopy or a single continuously tracked pole from 7/8 to 8/8.

Gate D: PARTIAL. Composition effects and modal tracking are supported; the physical homotopy portion is BLOCKED by discrete model semantics.

## 6. Gate E - exact network closure

The installed repository/API does not expose exact K, D, and Q objects with documented dimensions, units, construction path, determinant identity, or Q -> -1 boundary semantics. The audit therefore marks all five requested closure items BLOCKED. No graph Laplacian, arbitrary proxy, or fabricated identity was substituted.

Gate E: BLOCKED.

## 7. Claims and interpretation

The defensible result is narrower than a universal mechanism claim:

1. The original frozen discovery set contains an exact 255+1 robustness blocker pattern under its stated conditions.
2. Four high-PLL nonlinear TDS cases support the ordering V8 fragile, worst 7/8 intermediate, best 7/8 robust.
3. Aggregate-matched comparisons show composition-sensitive alpha differences, so converted MW/MVA alone do not explain every comparison.
4. The holdout census rejects the statement that the same H0 or H0.05 identity persists across all 24 designed conditions.
5. Independent finite-difference alpha validation is numerically unresolved/failing under the frozen step protocol because of strong eigenvalue sensitivity; this prevents a clean numerical-truth headline.
6. Exact closure theory and a physical continuous SG-to-GFL homotopy remain unavailable in this installed model.

The work that was intentionally not run is weak-node analysis, weak-link analysis, the 128-direction structured-radius campaign, repair optimization, planners, co-design, IEEE-68, EMT, and new controller search.

## 8. Final decision

Numerical audit: FAIL on the preregistered independent-alpha gate.
High-PLL TDS: PASS.
H0 over 24 holdouts: 42 total; counts [4, 2, 2, 0, 3, 0, 2, 2, 0, 0, 0, 19, 0, 0, 0, 1, 1, 2, 0, 3, 0, 0, 0, 1].
H0.05 over 24 holdouts: 47 total; counts [4, 2, 2, 0, 5, 0, 2, 3, 0, 0, 0, 19, 0, 0, 0, 2, 1, 2, 0, 4, 0, 0, 0, 1].
Penetration/composition/mixed verdict: composition-sensitive effects are present in aggregate-matched comparisons; no mixed intervention was run or claimed.
Mode mechanism: physical electromechanical/control modal distinction is supported, but continuous 7/8 -> 8/8 homotopy is blocked and the independent FD alpha audit is unresolved.
Closure: BLOCKED.
IAS poster redesign: NO.

FINAL CASE: C
