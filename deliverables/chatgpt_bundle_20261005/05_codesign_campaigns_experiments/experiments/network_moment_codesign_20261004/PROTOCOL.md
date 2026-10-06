# Exact network frequency moments and collective latency compensation

## Theoretical question

Can the cubic delay dependence of a locked type-II PLL be propagated through
the full grid's simple rotational singularity, and used as an exact constraint
on finite controller retuning? The all-PLL hidden block, not the earlier
singular speed-port hidden block, is used. No physical Laplacian is substituted
for the actual lossy network. Time-integral interpretations require stability.

The conjectured laws for arbitrary input/output maps are:
1. DC bus-frequency transfer is independent of PLL gains and delays.
2. At fixed Ki, changing Kp and delays leaves the signed area of every stable
   bus-frequency step-response difference zero.
3. Its first time-weighted integral is kappa*f_infinity, with
   kappa = l*diag(delta_tau/Ki-delta_Kp/Ki^2)*r/(l*K1*r).
4. The network-weighted scalar kappa=0 permits cross-site compensation.
5. This constraint can conflict with modal stability. A projected minimum-cost
   correction will be compared with sitewise delay compensation, at fixed rho.

The tools (Laurent residues, moment identities and convex QP) are classical.
Novelty is conditional on a useful full-network design consequence, not on
renaming a pseudoinverse or moment. No maximum replacement claim is targeted.

## Frozen experiment

Base: the preceding corrected design (rho30=rho37=.885, other rho=.875;
gains from interaction_decision_20261004/designs/corrected.toml). Original
gain limits stay unchanged. Base delays: all40ms. Ki remains fixed in design.

Delay patterns, fixed exogenously before comparison:
- Uniform +0.1, +0.5, +1.0ms.
- Heterogeneous +1ms at bus30; +1ms at bus37.
- Heterogeneous deterministic ramp0...1ms on buses30...39.

Three methods at each identical rho/delay:
unchanged gains; sitewise deltaKp=Ki*delta_tau;
collective moment equality plus critical-root margin correction.
Target real part −.06/s; required margin −.05/s. Positive-gain bounds unchanged.
The collective update minimizes squared relative Kp change, preserving the
exact scalar moment condition, with one active modal inequality. A trust
region is used and all catalog modes are checked; no global optimum asserted.

Analytical validation:
- all-PLL hidden-block regularity, gauge rank and nonzero derivative;
- 39 bus outputs and unit impedance-load parameters at buses8,16,29;
- exact Laurent coefficients from bordered recurrence versus small-frequency
  evaluation at s=.02,.01,.005,.0025 (1/s), plus gauge/projection error;
- test single-site Ki changes +10% at buses30,32,33,39;
- distinguish structural symmetry from exact independent binary64 entries.

If useful, full numerical DDE contour validation and five frozen nonlinear
events for the collective uniform+1ms design are the preregistered headline.
Sitewise comparator gets the same spectrum test; do not simulate an unstable
comparator merely to produce a dramatic plot. No event set/limit changes.
Heterogeneous delay cases receive spectral/numerical-moment tests only unless
a separately recorded addendum authorizes a correct heterogeneous integrator.

Negative findings and unresolved cases remain in the tables. Passing moments
does not imply a better nadir, RoCoF or stability. No physics uncertainty claim.
