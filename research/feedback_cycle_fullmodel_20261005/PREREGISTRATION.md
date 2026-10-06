# PREREGISTRATION — feedback-cycle falsification in the full IEEE-39 SG/GFL model

Written 2026-10-05 BEFORE any reveal computation of this campaign. Frozen by sha256 (see `PREREGISTRATION.sha256`).
Later edits are forbidden; deviations go to `NEGATIVE_RESULTS.md` as dated amendments.

## 1. Provenance (frozen)
- Repository HEAD `8f21ec0b3be48f217e280f615bfd1c619df4d7ff`, branch `codex/collective-interaction-bounds-20261003`.
- Tracked files dirty at start: README.md, THEORY.tex, poster main.tex, bnd_h4_mechanism files, `src/pd39/model.jl`. None is used by this campaign. Untracked files are many; only the inputs hashed in `derived/ENV_RECORD.json` are used.
- Julia 1.11.9; PowerDynamics 5.0.0; NetworkDynamics 1.3.0; OrdinaryDiffEqRosenbrock 2.7.1; SciMLBase 3.53.2 (from Manifest.toml). Python 3.11.10, numpy 2.0.1, scipy 1.14.0, pandas 2.2.2.

## 2. Canonical model (frozen, not modified)
- `experiments/nonlinear_codesign_20261001/ReducedDAE.jl` with `dc_convention=:physical_supply`, context `N.design_context(ROOT)`; exact linear export in `experiments/graph_gsp_codesign_20261003/model/*.csv` (Adev, Bv, Cs, Cf, Ds, Y, Etheta, Hv, Bp, Bi) and Python reader `experiments/interaction_decision_20261004/model.py` (`Model`, exact exponential delay, no Pade). Delayed nonlinear solver `experiments/graph_gsp_codesign_20261003/DelayedEvents.jl` (method of steps, Rodas5P).
- 204 differential states (SG, governor, AVR, converter current, DC and PLL states retained; network algebra eliminated exactly). SG/GFL sites: buses 30–39 (ten sites, order 30..39). Gauge: one rotational zero root, excluded.
- Operating point: initialized repository dispatch; P_GFL = rho·P0, P_SG = (1−rho)·P0 at each site.
- Baseline (frozen): rho_i = 0.875, Kp_i = 28.274333882308138, Ki_i = 246.74011002723395, uniform pure PLL error delay tau = 0.040 s.
- Stress case (frozen): same gains, uniform tau = 0.044 s ("S44", no retune). Historical result to be re-derived, not assumed: 8–10 roots beyond the margin; nonlinear events abort.
- All-SG reference: rho_i = 0 at every site (same topology, dispatch, loads). Modes that live only in decoupled GFL device states at rho = 0 are labelled GHOST and excluded from the all-SG physical spectrum.
- Gain bounds (frozen): 0.25–4 times nominal, nominal Kp0 = 2π·5, Ki0 = (2π·5)²/4, i.e. Kp in [7.8540, 125.6637], Ki in [61.685, 986.960]. rho is NOT a design variable here.
- Margin: sigma_req = 0.05 s⁻¹ (all non-gauge roots must have Re < −0.05). N_unstable = #(Re>0), N_margin = #(Re>−0.05), both non-gauge; never conflated.
- Frequency scan (linear return analysis): omega = 2π f, f in [0.05, 20] Hz, 4000 log-spaced points, along s = −sigma + jω with sigma = 0 and sigma = 0.05; refined around each located root.
- Candidate buses: 30..39. Preregistered sparse support S0 = {30, 33, 36, 37}.
- Frozen nonlinear events: the five admittance-load events of `DelayedEvents.CASES` = (bus, MW) in {(8,−100), (16,+100), (16,−100), (29,+100), (29,−100)}, 60 s horizon, pass iff |Δf| ≤ 0.5 Hz, RoCoF ≤ 0.5 Hz/s (0.5 s windows), V in [0.9, 1.1], SG governor/AVR slack ≥ 0.002.

## 3. Frozen external predictions (Python reduced model; NOT to be matched numerically)
P1 all-SG has no PLL/controller family comparable to the SG→GFL one. P2 SG→GFL replacement creates a distinct converter/synchronisation family (reduced model 7–8 Hz, frequency NOT required). P3 dangerous interactions are local to 35–36, 30–37, 33–34 (ranking is a prediction, not a constraint). P4 Kp-only has materially less useful authority than Ki or joint. P5 test S0 = {30,33,36,37} before any optimiser. P6 forcing near the limiting frequency excites PLL states more than 0.7f* and 1.3f*. P7 the phenomenon is distinguishable from an all-SG electromechanical resonance.

## 4. Operational definitions and falsification criteria (frozen)
- PLL-family mode: oscillatory mode (Im>0, f in [1,20] Hz) whose PLL-state participation fraction (sum of |l_k r_k| over PLL states theta, omega, xi of the ten GFL, normalised over all states) is ≥ 0.5.
- P1/P2: PASS if S44/baseline contain ≥ 1 PLL-family mode and the all-SG physical spectrum contains none; REFUTED if the all-SG spectrum has an oscillatory physical mode in [3,20] Hz with SG-state participation ≥ 0.5 within 0.5 Hz of a baseline PLL-family mode, or if no PLL-family mode exists at rho = 0.875. Branches are tracked by modal assurance criterion (MAC ≥ 0.7 step to step), never by frequency alone.
- Collective closure (return analysis): local factor regular if sigma_min(local factor) ≥ 0.1·median over the scan; collective closure if min over omega of sigma_min(I+Q) ≤ 0.05 while all local factors regular. Sign convention of Q to be verified by determinant parity at ≥ 20 random regular points, tolerance 1e-8 relative.
- P3: the dominant cycle is the pair with the largest |log|1 − p_ij|| contribution at the limiting root. P3 is PASS only if the top-2 pairs include at least one of {35–36, 30–37, 33–34}; PARTIAL if one of the three appears in the top 5 only; REFUTED otherwise. Whatever core dominates is reported.
- P4: PASS if the best Kp-only critical-root real-part shift per unit normalised gain effort (log-gain) is ≤ 0.5 times the Ki-only value, over the ten sites; REFUTED otherwise.
- P5/S0 repair: trust-region local LP/QP with exact-root recomputation after each finite step; step cap |Δ log K_i| ≤ 0.1 per iteration, ≤ 60 iterations; success iff N_margin = 0 (full-spectrum argument-principle count AND catalogue). Variants Ki-only, Kp-only, joint. Result reported whether positive or negative; no support change before S0 is finished.
- P6: PASS if peak PLL detector error at f* exceeds both 0.7f* and 1.3f* responses by a factor ≥ 3 at identical forcing amplitude on the stress case, and the repaired case suppresses the f* response by a factor ≥ 3 versus the stress case.
- Forced response (frozen before reveal): forcing = sinusoidal admittance-load modulation, peak 5 MW-equivalent (5 % of the event size), at the event bus among {8,16,29} with the largest linear input residue norm onto the limiting root (rule fixed here; computed once from the linear model), horizon 30 s, same solver and metrics as the events; frequencies 0.7f*, 1.0f*, 1.3f* where f* is the frequency of the limiting stress-case root at 44 ms. Cases: all-SG, S44, repaired S44.
- P7: PASS if the all-SG forced response at f* is not larger than that of the repaired case by a factor ≥ 3 and the S44 PLL-state growth is not explained by an SG electromechanical mode (PLL-state participation of the limiting root ≥ 0.5).
- RHP zeros: channels = each of the ten PLL detector outputs from each of the ten PLL injection inputs, plus the three event-bus load inputs to bus frequency; transmission zeros of the (delay-free-at-zero) square subsystems computed from the eliminated descriptor; reported only if confirmed with two independent methods. A negative result is recorded.
- Decision gate: STRONG_PASS requires distinct family (P1/P2) AND local factors regular with collective closure AND localised core (P3 PASS or PARTIAL) AND S0 or documented sparse support restoring N_margin = 0 AND forced-response/nonlinear agreement (P6 PASS and all five events pass). PARTIAL_PASS if some hold. REFUTED if the distinct family, the closure, or every sparse repair attempt fails with the above criteria. UNRESOLVED otherwise.

## 5. Not allowed
Changing the model, thresholds, events, delays, bounds, margin or candidate set after reveal; Pade in final root results; tuning to 7.879 Hz; forcing 35–36 as the winner; calling floating-point root counts "certified".
