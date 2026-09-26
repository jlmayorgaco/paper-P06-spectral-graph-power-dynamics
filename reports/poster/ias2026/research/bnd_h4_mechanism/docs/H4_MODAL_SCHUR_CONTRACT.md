# H4_MODAL_SCHUR_V1 contract

This is a narrow follow-up to the frozen Python TX4/H4 run. It does not
retune the model, alter P4, add graph modes, or overwrite the canonical run.

## Frozen input

- Baseline: `results/20260925T204017_549c3c07_h4_crossmode_v3`.
- Case: `H4=(30,33,35,37)`, `P4=(g,k,t,h)=(0.03625,1.425,1.5,1)`.
- Critical full-DAE pole: the positive-frequency transverse mode selected by
  the existing `mode_of` path.
- Retained operator: `F(s)=I+M(s)` from the exact retained port action space.
  The unnormalized `I+M` form is required so the eta=0 ablation retains
  physical local/self blocks; the normalized `I+Q` form would make the local
  diagonal identity by construction.

## M1 scalar Schur gate

At the full-DAE pole `lambda_c`, compute left/right singular vectors of
`F(lambda_c)` and use them as the first vectors of bi-orthogonal unitary
bases. In that basis partition `H=U^H F V` as a 1-by-1 critical block and a
7-by-7 complement:

```text
h_k(s) = H_kk(s) - H_kb(s) H_bb(s)^-1 H_bk(s)
```

Strict gate: `abs(h_k(lambda_c)) < 1e-8`, finite/invertible complement, and an
independent complex root of `h_k(s)=0` within `1e-6` of the full-DAE pole.
If the strict residual fails, the run stops before promoting M2/M3 claims.

## M2 collective-feedback ablation

Only after M1 passes, continue the same scalar root over
`eta=0,0.05,...,1` using

```text
h_k(s;eta) = H_kk(s) - eta H_kb(s) H_bb(s)^-1 H_bk(s).
```

The synthetic block operator scales both cross blocks by `sqrt(eta)` for
mode-tracking diagnostics. Store the root, real/imaginary parts, residual,
smallest singular value, and adjacent-root MAC at every eta. The primary
metric is `Re(lambda(eta))`; no damping interpretation is allowed.

## M3 derivative gate

Use

```text
d lambda/d eta = (H_kb H_bb^-1 H_bk) / (partial_s h_k)
```

and compare it with centered finite differences of the continued roots. The
derivative table is only claimable if the continuation remains valid and the
finite-difference comparison is finite and small.

## Outputs

All files are written only under one immutable run root. If M1 passes, the
only figures created are `F01_cross_feedback_root_locus`,
`F02_cross_feedback_complex_plane`, and
`F03_local_collective_boundary_overlay`. No graph modes, Taylor loops,
`p_min`, EMT, or new controllers are part of this run.

The exact Julia GFL11 reproduction is a separate parity check using
`code/tx4/run_tx4_exact_p4_ieee39.jl`; the historical 82-state GFL10 runner is
not used as a gate.
