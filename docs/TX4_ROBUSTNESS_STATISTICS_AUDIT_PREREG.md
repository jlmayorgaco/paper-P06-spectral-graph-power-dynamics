# TX4 Robustness Statistics Audit — Preregistration

Date: 2026-09-20  
Branch: `research/tx4-robustness-statistics-audit`  
Parent: `a9467a457fb029715945a149957af0e7dc9c6b23`  
Frozen exact parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`

## Immutable inputs

- V4 = {30, 33, 35, 37}; H4 is the full portfolio; H0 is the minimal-blocker
  antichain over all 16 subsets.
- Existing QMC and MC coordinates, seeds, uncertainty bounds, model,
  equilibrium solver, central-difference linearization, and 0.3–1.5 Hz EM
  mode band remain fixed.
- No final paper/poster PDF is an input to the corrected Vancouver package.
- No push is authorized. The original robustness branch and frozen exact
  worktree are out of scope for modification.

## Definitions

For a condition c, let `alpha_EM(S,c)` be the maximum real part of the
non-gauge eigenvalues in the registered 0.3–1.5 Hz band for portfolio S.

`U0(c) = {S : alpha_EM(S,c) >= 0}`.

`H0(c)` is the inclusion-minimal members of U0: a blocker H is retained only
when no strict subset of H belongs to U0. The primary endpoints are:

- `H4_UNSTABLE`: H4 belongs to U0;
- `H4_PRESENT`: H4 belongs to H0;
- `EXACT_H4`: H0 contains exactly one blocker and it is H4;
- `NONCOMPOSABLE`: some H in H0 has cardinality at least two;
- `KAPPA`: minimum blocker cardinality, or NULL when H0 is empty.

The following implications are asserted and tested: EXACT_H4 implies
H4_PRESENT implies H4_UNSTABLE, and EXACT_H4 implies NONCOMPOSABLE.

The prior source code contains a primary threshold of 1e-8 s^-1 for legacy
flags, but no registered U005/H005 threshold. U005/H005 will therefore be
reported only if an exact pre-existing engineering threshold is located in
the audit trace; otherwise they are explicitly marked not analyzed.

## Statistical handling

QMC results are coverage summaries over the fixed 4096-point scrambled Sobol
design and are not called probabilities. MC results use the fixed 5000-point
PCG64 independent-uniform design and may be described as assumed bounded-box
probabilities with Wilson intervals. KAPPA is summarized as a PMF. H4_UNSTABLE,
H4_PRESENT, and EXACT_H4 are separate endpoints.

If the surrogate gate fails, the exact QMC/MC fallback is mandatory: 4096 and
5000 existing coordinates, all 16 portfolios, 145,536 portfolio evaluations,
checkpointed every 250 conditions, with no sample regeneration or reduction.

## Required provenance

Every master row is assigned one of EXACT_DAE,
SURROGATE_EXTRATREES, DERIVED_FROM_EXACT, or UNKNOWN. Surrogate-derived
minimality is prohibited from primary reporting unless all strict gates pass.

