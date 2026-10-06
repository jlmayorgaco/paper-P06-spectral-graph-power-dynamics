# Claim ledger — supported statements only

## Claims supported by completed calculations

- **C0 — No-delay implementation parity — NUMERICALLY_VALIDATED.** The 87.5% GFL seed, 90.0470% historical joint candidate, and 89.8604% fixed-gain candidate reproduce the no-delay trim, 203-pole spectrum, rightmost pole, and tested no-delay event outputs between ReducedDAE and independently compiled PowerDynamics within the tolerances in `baseline_reproduction/D00_GATE_REPORT.md`. This is parity, not a security-feasibility claim.
- **C0a — Existing joint design misses two external frequency limits — NEGATIVE_RESULT.** The reproduced 90.0470% candidate reaches 0.50569 Hz for bus 8, −100 MW and 0.51772 Hz for bus 16, −100 MW against the frozen 0.5-Hz limit. It therefore is not a fully feasible lower-bound point, even at zero added delay.
- **C1 — Exact retarded characteristic form — EXACT_IDENTITY for the declared fixed-support delayed-error model.** The implemented matrix is Δ(s)=sI-A0-Σ Ai exp(-sτi), with the PLL detector error delayed only in the Kp/Ki branches. Its τ=0 spectrum matches the finite ODE spectrum. This mathematical form alone is not a computed complete DDE spectrum.
- **C2 — Local root sensitivities — NUMERICALLY_VALIDATED / SUPPORTED_LOCAL.** The tracked dominant zero-delay conjugate pair's τ, Kp and Ki derivatives pass 27 centered finite-difference checks at 3 frozen patterns and 3 ports. This does not include ρ and does not establish full DDE root completeness.
- **C3 — Lossless branch Taylor expansion — NUMERICALLY_VALIDATED** only for the explicitly declared sinusoidal lossless branch-power map and the 60 test directions/amplitudes in `TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv`. At the largest tested amplitude, 0.1 rad, median relative error falls from 4.37e-4 at first order to 2.28e-4 at second order and 5.30e-8 at third order. This is not a lossy AC or nonlinear stability result.
- **C4 — Delay vectors fixed before evaluation — NUMERICALLY_VALIDATED.** The same ten delay values (mean 25 ms, standard deviation 16.8203 ms, maximum 50 ms) were frozen into 104 bus assignments and hashed before inspecting delayed roots. This is experimental design provenance, not evidence that assignment changes maximum replacement.
- **C5 — Exact Schur damping operator — BLOCKED.** After gauge deflation, the proposed hidden block is singular at s=0 (condition estimate ~3.10e18), so the retained Schur operator and its slope are undefined for this partition. The non-deflated table is preserved only as a discarded diagnostic. No D_G, L_D, or L_tau physical claim is supported.

## Questions that remain unresolved

- **Q1 decomposition (D_G=D_{base}+L_D-L_\tau): BLOCKED.** The required gauge-deflated Schur slope is undefined for the current partition, so none of the three physical terms is accepted. The algebraic decomposition remains only a proposed definition.
- **Q2 first-order delay approximation: INCONCLUSIVE.** The simplified PLL approximation has not been compared against a full-model critical damping loss over a validated delay interval.
- **Q3 same multiset, different placement changes maximum replacement: INCONCLUSIVE.** No complete DDE spectrum or placement-specific co-design was obtained.
- **Q4 commutator predicts replacement: INCONCLUSIVE.** χτ and roughness are computed from the pre-frozen graph, but no validated replacement capacity exists to correlate.
- **Q5 GSP gain compression: INCONCLUSIVE.** No graph-mode gain-restricted optimization was run.
- **Q6 reduced analytical solution predicts active SG retention: INCONCLUSIVE.** No delay-specific full design was solved.
- **Q7 full optimum has at most two partial SGs: NOT TESTED.** The two-anchor result is only a theorem for the specified two-constraint reduced LP; it does not imply the full event-constrained design structure.
- **Q8–Q10 maximum replacement, optimal gains, and delay shadow prices: NOT ESTABLISHED.**
- **Q11 (R_F\le R^\star\le R_U): BLOCKED.** There is no validated feasible delayed design or rigorous upper bound.

## Poster decision

**POSTER READY: NO** for a claim about a delay-dressed maximum SG→GFL replacement frontier. The blocking gate is `BLOCKED_EXACT_DDE_SPECTRUM`; nonlinear delayed events, delayed co-design, and the upper-bound calculation are consequently absent. The current evidence supports a methods-progress poster only if clearly framed as a failed/blocked validation study, not as a solved replacement frontier.
