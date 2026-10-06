# Physical collective damping: preregistered mechanism test

Frozen before executing the new calculations, 2026-10-03. This is a mechanism
experiment, not a new replacement optimization or a full event-security campaign.
Existing frozen outputs are read-only. No commits or pushes.

## Question and falsification

Does the latency-sensitive PLL mode have a material, physically interpretable
collective damping contribution at the retained synchronous torque/speed ports?
If so, does a local nonlinear delayed perturbation reproduce its predicted decay
or growth? Merely obtaining a dense matrix, a Schur identity, or an improvement
with retuning is not evidence of a novel method or greater replacement capacity.

## Frozen inputs

- Three existing designs: N_nominal, Z_zero_delay_tuned, full20_step15_medium.
- Their existing rho values and initialized dispatch, without optimization.
- Uniform PLL error delays 0, 20, 30, 40, 45 ms. Additional evaluations only for
  continuation and locating Re(s)=0 crossings of the same tracked roots.
- Frequency scan: 120 log-spaced values 0.01--20 Hz, plus tracked mode frequencies.
- Full state dynamics with physical_supply DC convention. No Pade approximation.
- Keep SG speed perturbations at buses 30--39 as physical port outputs, with
  externally added mechanical torque expressed on retained rating as inputs.

## Gates and metrics

1. Equilibrium infinity residual <1e-7; physical Jacobian agrees with independent
   full-RHS automatic differentiation to 1e-8 relative; eigenvalue sets agree
   with existing characteristic (after gauge removal) to 1e-5 absolute.
2. Finite-frequency hidden block nonsingular. Report conditioning. Check port
   Schur/resolvent equality and physical harmonic-work identity to 1e-7 relative.
3. Exact exponential characteristic root residual <1e-8. Roots are tracked local
   branches, not an exhaustive certificate for the infinite DDE spectrum.
4. Report full Hermitian port damping, its nodal diagonal and its off-diagonal
   contribution on the actual mode. Materiality: deleting the off-diagonal term
   changes the sign of mode-projected harmonic damping OR its magnitude exceeds
   10% of the sum of absolute diagonal and off-diagonal contributions.
   This deletion is an attribution diagnostic, not a realizable controller.
5. Nonlinear method-of-steps validation for N and full20 at 40 and 45 ms, plus
   zero delay. Initialize a tracked eigenmode history with max PLL-angle
   amplitude 1e-5 and 1e-4 rad; simulate 3 s. Fit mode growth using the matching
   complex left eigenvector. Compare nonlinear and linear DDE histories/outputs.
   Numerical convergence: repeat at tighter tolerance and shorter maximum step
   if baseline disagreement exceeds 2%; report amplitude effects separately.
   This is local nonlinear validation, NOT the five frozen +/-100 MW events.

## What cannot be claimed

No new maximum penetration, global optimum, all-event security, zero-frequency
Laplacian law, or novelty proof follows from these tests. A negative Hermitian
impedance direction is a harmonic incremental-work diagnostic; it does not alone
prove closed-loop instability. An input/output elimination can hide modes when
its eliminated block is singular. No passivity theorem without additional gates.

## Evidence statuses

EXACT_IDENTITY, PROVED_REDUCED_MODEL, NUMERICALLY_VALIDATED, SUPPORTED_LOCAL,
SUPPORTED_EMPIRICAL, NEGATIVE_RESULT, BLOCKED. Statuses must remain explicit.
