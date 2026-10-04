# Adaptive follow-up, frozen before evaluating compensated replacements

The original fixed grids produced no valid false acceptance and five false
rejections of individual modal-margin predictions. 112 of 550 comparisons were
invalid after tracking collisions (mostly common Kp). Preserve all these results.
The intended exclusion story has therefore NOT been obtained by that campaign.

New hypothesis: individually pole-preserving SG replacement/PLL retuning may
fail to compose when several sites implement their changes simultaneously.
This tests an engineering action, not arbitrary gain damage. The selection is
explicitly adaptive discovery, not a held-out confirmation.

Anchor: family I, k=7 from the frozen sweep: rho=.875, original Kp, Ki*1.14,
tau=40ms at every site. Its leading PLL mode is still left of -.05/s; choose
that mode (historical mode 2) as the target. This is a numerical anchor, not a
nonlinear-secure or globally optimized design. Other roots still require checks.

For each site and delta rho in {.0025,.005,.01,.02,.04}, change ONLY that rho.
At the fixed target lambda, open only that PLL channel in the exact full DDE:
D_-i=lambda I-A0-sum_(j != i) Bj exp(-lambda tau_j) Cj.
ap=exp(-lambda tau_i) Ci D_-i^-1 Bp_i and similarly ai with Bi_i.
Solve [[Re ap, Re ai],[Im ap, Im ai]] [Kp_i,Ki_i]^T=[1,0]^T.
Retain only solutions within the historical physical gain bounds and verify
the full characteristic residual. This is exact single-pole compensation,
not an optimum nor preservation of the rest of the spectrum.

Combine all 45 site pairs, plus the all-ten combination, at each increment.
Evaluate every valid combination, retaining failures. For selection require a
joint root violation >=1e-4/s while the matched singleton target retains its
original margin. Track through at least 8 intermediate replacement steps and
record pole moves; collisions require a cluster treatment rather than relabeling.
Prioritize the smallest replacement increment with a qualifying pair, then the
largest violating real part, breaking ties by ascending bus pair.

Attempt a finite, continuous proof for that pair using the exact interaction
determinant, singleton root counts, signed two-step correction and fourth-order
pair-walk remainder. If inconclusive, retain the result as numerical only.
If every pair is harmless, this follow-up is negative; do not redefine success.
All-ten changes are secondary; no assertion about physical transmission cycles.
