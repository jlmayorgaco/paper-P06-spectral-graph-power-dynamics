# F0 — Architecture-invariant voltage-frequency metric

**Gate: PASS_CONTINUOUS.** The audit finds that ExpG used an inertia-weighted retained-SG COI output and ExpP used retained-SG COI at the retained machine; neither freezes an architecture-invariant bus-frequency sensor or filter bandwidth. ExpQ's Savitzky–Golay bus derivative was a diagnostic, not a frozen primary requirement.

For design measurements, this stage uses the same passive quantity for every architecture: unwrapped busbar-voltage phase at generator buses 30–39, differentiated at a fixed 100 Hz sample rate. The 1 MW bus-16 event is evaluated for all ten buses at `rho=0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4`; all-SG is the reference. A sampled unfiltered finite difference is retained along with cubic Savitzky–Golay phase differentiators with 0.04, 0.10, 0.20 and 0.40 s supports. The filter outputs are a bandwidth sensitivity study, not a post-hoc choice of the passing metric.

At `rho=0`, the algebraic bus-voltage phase has a jump of up to 7.536547954173718e-6 rad for the 1 MW step. Its ideal unfiltered derivative contains an impulse; consequently the continuous-time unfiltered RoCoF supremum is unbounded. Sampled and filtered RoCoF values are finite but depend on estimator bandwidth. The inherited 0.5 Hz/s design limit was previously applied to SG COI and has no frozen bus-measurement bandwidth; its application to bus-frequency RoCoF is therefore **METRIC_DEPENDENT**.

Continuity uses condition-aware tolerance `max(1e-8 × scale, 100 eps × max(cond(A),cond(A_all-SG)) × scale)`. The rho=1e-8 endpoint maximum absolute discrepancy is 5.618558049369504e-11 versus tolerance 2.9970065065933198e-5; pass is true. The metric is continuous at the tested GFL-insertion boundary only to this finite-horizon, sampled linear test. See the detailed tables for each output/window, rho level, and conditioning.

No PowerDynamics call or optimization was used. Evaluation cost: 1000 per-output continuity comparisons over 60 architectures; elapsed 414.4370000362396 s. This gate allows proceeding only if the endpoint continuity checks pass. It does not resolve a defensible primary RoCoF bandwidth.
