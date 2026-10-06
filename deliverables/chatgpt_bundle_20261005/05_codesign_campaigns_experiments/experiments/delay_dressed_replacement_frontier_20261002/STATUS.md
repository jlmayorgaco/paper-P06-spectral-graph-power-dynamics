# Campaign status — 2026-10-02

**Overall: BLOCKED for the requested delay-dependent maximum replacement frontier.**

The no-delay baseline is numerically reproduced independently in ReducedDAE and PowerDynamics. The exact retarded characteristic matrix and its τ=0 identity are implemented, but a complete rightmost-root count for positive delay was not established. Consequently, no delayed candidate was optimized or called feasible, and no value of (P_{\mathrm{GFL}}^{\max}(\tau)), (\rho^\star), (K_p^\star), or (K_i^\star) is reported.

Completed evidence:

- **D0 — NUMERICALLY_VALIDATED:** baseline trim/equilibrium, 203 finite poles, rightmost pole, and design plus five frozen no-delay events agree across ReducedDAE and independently compiled PowerDynamics within declared tolerances. See `baseline_reproduction/D00_GATE_REPORT.md` and `baseline_reproduction/TABLE_D00_*`.
- **No-delay feasibility warning:** parity is not security feasibility. The reproduced 90.0470% historical joint candidate exceeds the 0.5-Hz frequency limit in two of five holdouts: 0.50569 Hz at bus 8, −100 MW and 0.51772 Hz at bus 16, −100 MW. Its measured RoCoF and voltage metrics pass; actuator slack was not recomputed in this D00 event-parity run. It cannot be counted as a fully feasible replacement lower bound.
- **D1 — BLOCKED_EXACT_DDE_SPECTRUM:** exact characteristic matrix assembled without Padé; τ=0 spectrum identity passes; 3 delay patterns have a locally tracked dominant ODE conjugate pair with small residuals. The rightmost DDE root count remains unresolved, so these are not a complete DDE spectrum.
- **D2 — SUPPORTED_LOCAL / overall BLOCKED_PARAMETER_SET_INCOMPLETE:** τ, (K_p), and (K_i) simple-root derivatives pass 27 centered finite-difference checks for the tracked pair. ρ sensitivity and complete-root-set validation remain outstanding.
- **D3 — NUMERICALLY_VALIDATED:** lossless sinusoidal branch-power Taylor tensors pass 60 deterministic direction/amplitude checks. This is not a lossy AC or nonlinear stability result.
- **D4 — BLOCKED_SCHUR_SINGULAR_AFTER_GAUGE_DEFLATION:** after explicitly projecting out global rotation, the hidden block Δ_cc(0) is singular for the declared retained SG-angle coordinates. The Schur slope is undefined under this partition; the earlier no-gauge calculation is preserved only as a discarded diagnostic and is not evidence.
- **Delay placement:** all 104 same-multiset delay vectors were frozen and hashed before delayed roots were inspected. The L_P used for graph descriptors is the valid lossless graph Kron reduction; the exact lossy angle Jacobian is nonsymmetric and is not called a Laplacian. No descriptor-to-capacity relation was evaluated.

Not completed because the D1 and D4 gates block full-design validation:

- delayed full-spectrum frontier/co-design and comparisons against common/fixed/nodal gain strategies;
- full graph-mode gain compression and placement-capacity correlations;
- nonlinear method-of-steps validation of delayed events;
- delay shadow prices and a validated upper bound / optimality gap.

No Padé poles or trajectories are used as primary truth. Existing frozen files were not changed. No commit or push was made.
