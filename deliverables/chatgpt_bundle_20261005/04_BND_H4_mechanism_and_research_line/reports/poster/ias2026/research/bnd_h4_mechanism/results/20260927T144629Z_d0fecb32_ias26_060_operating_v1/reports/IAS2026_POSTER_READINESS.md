# IAS2026 poster evidence readiness — 20260927T144629Z_d0fecb32_ias26_060_operating_v1

Decision: **READY_WITH_EXPLICIT_LIMITATIONS**

Monte Carlo: 967/1000 valid; 33 physically infeasible; zero replacements.
TDS: 60/90 completed with a finite primary fit; 56/60 primary-estimator sign agreement; 4 disagreements, all 4 at 0 Hz with 4 unavailable R2 diagnostics, so mode identity is unresolved. 30/90 solver failures are retained.
Cross-code: same-model Python–Julia reproduction; mode identity unresolved.

## Safe claims
- Nominal H4 is an inclusion-minimal transverse blocker: all 15 proper subsets stable; H4 unstable.
- At the frozen nominal boundary, physical local factors remain nonsingular while collective closure approaches singularity.
- H4 minimal-blocker pattern occurs in 190/967 valid frozen synthetic scenarios.
- Fixed g=0.25 improves α⊥ in 594/967 and rescues H4 in 186/967 valid scenarios.
- Python–Julia verdict agreement 100/100; the TDS primary estimator agrees in sign in 56/60 completed finite fits.

## Claims not supported
- MC event frequencies are real-world IEEE-39 probabilities.
- The H4 blocker persists for every perturbed operating point or every grid.
- g=0.25 universally stabilizes H4 or never deteriorates stability margins.
- Mode identity is tracked across scenarios or Python and Julia eigenmodes are physically independent validation.
- The nonlinear phasor DAE trajectories are EMT/field validation.
- IAS26-010 M1 strict gate has passed or any eta/M2/M3/GFM claim is established.

F1 and F2 are exact copies of audited frozen assets. P3 is the frozen MC figure. P4 combines paired fixed-intervention MC results with only the three completed canonical D2/A1 phasor-DAE traces. The ensemble is synthetic, not a real-world probability sample.
