# F0–F7 experiment program — gate outcome

**PROGRAM_STATUS: PARTIAL — stopped at the pre-registered F2 modal gate.**

The two pasted requests were identical. Experiments A–E and their existing
artifacts were left untouched. F0 and F1 passed. F2 showed that the physical
conductance graph does not diagonalize the actual detailed GFL/network pencil;
the spectral gain-law and graph-controller branches were therefore stopped.
F2 also retained exact coupled-mode Schur self-energy diagnostics for a
separate pathway-design direction. The analytic stages F0–F2 did not import
or call PowerDynamics.

## Gate evidence

- **F0 PASS:** The pre-registered graph is the Hermitian part of the complex
  passive branch Ybus after Kron reduction to generator ports 30–39.
  (L_c=\operatorname{Re}(Y_{port})) is SPD with eigenvalues 0.0011275 to
  5.01334 p.u. conductance, condition number 4446.27, and no artificial
  grounding. Direct and Kron port currents agree to relative residual
  (4.08\times10^{-16}). The reduced susceptance is indefinite; ExpC's C33
  synchronizing backbone also has a negative eigenvalue and is not square-rooted.
- **F1 PASS:** For each isolated commuting mode, (d_k^*=2\sqrt{\nu_k}).
  A deterministic five-mode SPD realization gave maximum analytic-to-numeric
  pole error (2.79\times10^{-7}\,\mathrm{s}^{-1}) and modal-transform
  residual (1.24\times10^{-14}). This is a modal reference result, not a
  global noncommuting-damping claim or a GFL result.
- **F2 NON_MODAL:** In the preregistered F0 graph basis, conductance coupling
  is (5.5\times10^{-17}), but the reactive off-diagonal ratio is 0.2914,
  above the 0.05 gate. The heterogeneous port transfer ratio is 0.2837 and
  the full descriptor residual is 0.0614. Hence the detailed model does not
  admit the proposed scalar (F(s;\nu,K_p,K_i)=0) reduction in this basis.
- The exact SimpleGFLDC PLL limit retains a cubic
  (\tau s^3+s^2+VK_p s+VK_i=0); the classical quadratic follows only as
  (\tau\to0). The isolated PLL limit itself contains no graph eigenvalue.
- F2 computed per-mode exact Schur Γ with zero reconstruction residual at
  the frozen C0 critical pole. It did not apply a gain correction because no
  graph coefficient seed exists.

## Not run after the stop gate

F3 modal gain optima and law fitting; F4 ideal/filter/node-local controllers;
F5 controller correction and recovery ratios; F6 graph-class replacement
capacity and KKT audits; and F7 hash-frozen PowerDynamics and time-domain
validation. Their reports say `BLOCKED` or `DIAGNOSTIC_ONLY`; no missing
performance cells are filled from unrelated ExpD/ExpE results. Current-limit
safety remains unavailable because SimpleGFLDC has no hard current limiter
or validated current rating.

## Required terminal summary

```text
GRAPH_PLL_PROGRAM_STATUS: PARTIAL
PHYSICAL_GRAPH: PASS
MODAL_STRUCTURE: NON_MODAL
GAIN_LAW_CLASS: NOT_APPLICABLE_F2_STOP
DISCOVERED_KP_LAW: NONE
DISCOVERED_KI_LAW: NONE
IDEAL_GRAPH_CONTROLLER_PARAMETERS: NOT_RUN
LOCAL_GRAPH_FILTER_ORDER: NONE
LOCAL_GRAPH_PARAMETER_COUNT: NOT_RUN
GRAPH_PERFORMANCE_RECOVERY: NOT_RUN
SELF_ENERGY_CORRECTION_GAIN: NOT_RUN
FREE_GAIN_REFERENCE_PERFORMANCE: NOT_RUN
UNIFORM_REFERENCE_PERFORMANCE: NOT_RUN
MAX_GFL_MW_UNIFORM: NOT_RUN
MAX_GFL_MW_GRAPH: NOT_RUN
MAX_GFL_MW_SE_GRAPH: NOT_RUN
MAX_GFL_MW_FREE: NOT_RUN
CAPACITY_RECOVERY_GRAPH: NOT_RUN
CAPACITY_RECOVERY_SE_GRAPH: NOT_RUN
POWERDYNAMICS_MAX_POLE_ERROR: NOT_APPLICABLE_NO_FROZEN_CANDIDATE
ROCOF_FINAL: NOT_RUN
CURRENT_HEADROOM_FINAL: NOT_AVAILABLE_HARD_LIMIT_NOT_MODELED
MAIN_THEORETICAL_RESULT: The isolated commuting second-order mode is maximally damped at d*=2*sqrt(nu).
MAIN_POWER_SYSTEM_RESULT: The passive conductance port graph is SPD, but reactive coupling prevents useful graph-modal separation of the detailed GFL network.
MAIN_LIMITATION: The 29.14% reactive cross-mode coupling exceeds the pre-registered 5% gate.
PUSH: NO
```
