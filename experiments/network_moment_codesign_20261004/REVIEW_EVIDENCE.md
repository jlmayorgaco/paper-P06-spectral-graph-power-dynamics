# Independent evidence review — 4 October 2026

Scope: read-only review of the exported maps, Laurent implementation, frozen
protocol/addenda, retained failed runs, and the running multimode co-design.
Only this review file was written. No trajectories or heavy simulations ran.

## Verdict

The moment recurrence, compensation sign, stable resolvent-difference formula,
and full-network event/output export are correct under the stated structural
assumptions. Independent finite-difference output checks and a full 204-state
transfer comparison passed. The results support an exact low-frequency
structural law with numerical validation; they do not establish a broadband
predictor, nonlinear performance law, or stability guarantee from the moment
constraint alone.

The co-design run was incomplete at this review snapshot: TABLE_06 contained
six rows, covering uniform +0.1 and +0.5 ms only, out of 14 frozen patterns.
Those collective candidates satisfy the eleven-root catalog margin, whereas
the sitewise candidates do not. Full-spectrum and nonlinear gates remain open.
This review must not be read as approval of subsequently generated results.

## Independent compact checks

The independent full-system frequency transfer was evaluated as

    H_full(s) = (1-exp(-s W))/(2 pi W) *
        [O Delta(s)^(-1) (Binput + B exp(-s tau) detector_input) + phase_input].

It was compared with Moments.transfer(..., restore=False), using all 39 bus
outputs and all three unit-MW load parameters. Relative discrepancies were:

| Positive real Laplace sample s, 1/s | Relative full/reduced discrepancy |
|---:|---:|
| 0.02 | 1.6593e-9 |
| 0.01 | 3.6832e-10 |
| 0.005 | 7.3782e-10 |
| 0.0025 | 1.2498e-8 |

These are positive real Laplace samples, not sinusoidal frequency-response
samples. The formulas include direct algebraic phase feedthrough.

Additional checks:

- Exported detector times physical gauge: norm 9.97e-15.
- Exported phase-state map times physical gauge minus ones: norm 2.10e-14.
- Raw A0 times physical gauge: norm 1.83e-9, before state/scale normalization.
- Hidden dynamics have exactly zero omega/xi columns. The discarded PLL rows
  of Binput have norm 4.08e-15, and phase output omega/xi columns are zero.
- Independent centered load finite differences of all-bus phase, at step
  0.001 MW, agreed with phase_input to relative errors 1.084e-9, 1.808e-9,
  and 1.228e-9 at buses 8, 16, and 29; maximum absolute error was 3.02e-13.
- Independent centered state finite differences, step 1e-6, agreed with the
  phase-state columns for Julia states 12 and 166 to relative errors 9.85e-10
  and 2.03e-9. State 15 is a structurally zero phase-output column; its finite
  difference was zero to 8.31e-13 absolute error (relative error is undefined).
- The existing source input-parity checks have maximum relative error 2.17e-9.
- Existing pole-gradient finite differences have maximum relative error
  3.99e-7 across the six reported Kp/delay directions.
- MOMENT_INPUT_LOCK_V1 and MULTIMODE_SOURCE_LOCK_V1 matched the inspected code
  and inputs at the review snapshot.

The full 78D event solve in export_inputs.jl correctly avoids losing the
event-dependent passive-bus lift contribution. Its model is the corrected
nonuniform-rho base, not the earlier uniform-rho model.

## Mathematical scope and numerical qualification

The coefficient recurrence is consistent with K(s) x(s)=g(s),
x(s)=x_minus1/s+x0+s*x1+..., and the 0.5-second bus-frequency window. Its
bordered solves enforce the required gauge solvability conditions. The
stable difference implementation correctly computes H_b-H_a without directly
subtracting nearly equal transfer matrices.

At fixed rho, positive Ki, gain-independent input/output maps, hidden-block
regularity, rank(K0)=9, and nonzero d=ell*K1*r:

- H0 is independent of Kp, Ki, and delay.
- At fixed Ki, changes in Kp/delay preserve H1.
- Delta H2 = -kappa H0, with the specified signed network weights.
- Therefore the difference of stable unit-step responses has zero signed
  area and first time-weighted integral +kappa H0.

The time-integral statements require stability and integrability. They cannot
be assigned to an unstable sitewise comparator merely because its formal
Laurent coefficients exist. The exported inputs describe the linearization
with respect to impedance-load parameters; exact identities for these
transfers are not exact identities for finite nonlinear +/-100-MW events.

Reported coefficient checks have maximum relative errors 3.60e-6 for the
cubic change and 1.85e-10 for the quadratic Ki change; lower coefficients are
unchanged in the implementation. These recurrence checks reuse the same
structural coefficients and are algebraic implementation checks. The full
transfer and independent output finite differences above supply additional
independent validation.

The finite-s asymptotic error at s=0.0025 is still about 6.616% for the uniform
delay cases and 6.972% for the ramp. Do not call this a broadband predictive
law or infer good approximation at the approximately 5-Hz oscillatory modes.

The code explicitly restores G0*r=0 and T0*r=ones. Corrections have norms
6.69e-12 and 1.29e-11. The raw versus restored difference is small in the tested
domain (roughly 6.1e-8 relative to the predicted coefficient at s=0.0025),
but this remains a symmetry-restored numerical germ, not an interval proof
for every independently rounded CSV entry. Raw hidden condition number is
3.89e9 and the balanced value remains 2.22e7. Preserve these disclosures.

The network weights are signed and coordinate-free moment sensitivities.
They should not be described as probabilities, positive participation shares,
or graph-Fourier modes without a separately demonstrated graph connection.

## Co-design comparison and outstanding gates

The inspected multimode QP uses the correct equality in relative-gain
coordinates, all eleven inherited catalog-root inequalities, physical gain
bounds, and a local trust region. Adaptive predictor/corrector continuation
rejects detected root collisions. Final candidate roots are reevaluated.
Its norm objective supports a local minimum-effort search description, not
a global minimum-cost theorem or exhaustive stability optimization.

Current snapshot:

| Pattern | Sitewise catalog alpha | Collective catalog alpha | Sitewise relative-gain norm | Collective relative-gain norm |
|---|---:|---:|---:|---:|
| Uniform +0.1 ms | -0.0033543 | -0.0600001 | 0.0031324 | 0.0051475 |
| Uniform +0.5 ms | +0.2196771 | -0.0600001 | 0.0156620 | 0.0382840 |

This supports the bounded observation that preserving the moment sitewise
does not preserve the catalog margin, while using the remaining controller
freedom can restore that margin in these cases. It does not show lower effort
than sitewise tuning. Both rules share the replacement, delay, Ki, and bounds;
the collective rule additionally uses a stability target. The comparison is
therefore a comparison of these specified rules, not universal superiority
over local control. An optional no-moment-constraint modal QP would quantify
the effort cost of preserving the moment, but is not required for the narrower
existence claim.

The single-mode failure is correctly retained: its +0.5-ms collective result
has an unstable inherited root near +0.1770/s. The addendum accurately labels
the revised multimode method and all-site placement check as exploratory.

Required before campaign-level claims:

1. Complete the frozen grid or save explicit unresolved/failure rows. Current
   main() calls collective() before either comparator; an uncaught collective
   failure aborts the rest of the grid and skips both comparator rows for that
   case. Failure must not be silently interpreted as a negative or omitted.
2. Retain catalog-only wording until independent complete-region DDE checks
   and relevant high-frequency coverage pass. Eleven roots are not a complete
   spectrum. A small characteristic residual is not a root-count certificate.
3. Run the registered nonlinear checks only for a valid headline candidate.
   Heterogeneous-delay cases require a correct heterogeneous integrator and
   the stated additional record before nonlinear claims.
4. Save a final provenance manifest binding moment code/exports, corrected
   base, imported old model.py and its source matrices, inherited root catalog,
   protocols/addenda, dependency versions, and final outputs. The present
   individual source locks do not bind this whole dependency chain.

No mutation of experiment code, data, designs, or earlier results was made
during this review.
