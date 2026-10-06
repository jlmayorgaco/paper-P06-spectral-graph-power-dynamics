# Claim status vocabulary and honest novelty assessment

Status labels used across the repo: EXACT_IDENTITY, EXACT_CONDITIONAL_CLOSED_FORM, PROVED_REDUCED_MODEL, GLOBAL_OPTIMUM_REDUCED_MODEL, LOCAL_FULL_MODEL_KKT, FULL_MODEL_VALIDATED, NUMERICALLY_VALIDATED, SUPPORTED_LOCAL, CERTIFIED_FULL_MODEL, NEGATIVE_RESULT, BLOCKED. A reduced-model optimum is never promoted to a full-model global optimum.

## Math re-check done on 2026-10-05 (independent, numeric)
- Pair determinant identity: verified on random matrices (14.312434077 on both sides) and by derivation (factor diag(B_i,B_j), Schur complement).
- Interaction contour formula: z + z0 - zi - zj equals the first log-derivative moment of the determinant ratio when root counts are equal.
- Herm(rank-one) has inertia (1,1) unless the two vectors are parallel (checked numerically).
- PI gain map: Kp*lambda + Ki = W split into real/imaginary parts — correct.
- Two-node floor: area law A_1 = P m/4 - asin(P/2)/2 matches the nonlinear simulation to 6 decimals; zero crossing at m = pi/3. Area >= 0 is necessary, not sufficient (m = 1.06 has positive area yet bus 1 overshoots by 0.148).
- Defect fixed: a poster panel said "no local retune can move" the floor; the theorem only covers inertia redistribution and zero-net-energy supplementary control. Notation clash: m = number of enclosed roots (interaction) vs total inertia (floor).

## Novelty (limited literature search, NOT proof of priority)
Known/classical: interaction measures in decentralized control (Rosenbrock, RGA, Gershgorin); multi-converter PLL interaction and grid-structure effects on PLL synchronisation (e.g. arXiv 1903.05489); PI pole placement with delay; inertia placement (Poolla-Bolognani-Dorfler); topology effects on frequency response; integral/area constraints on overshoot (Aalipour-Das); aggregate-neutral damping (Kohlhaas-Kotyczka 2026); Coletta-Jacquod PRE 2016; Gorbunov et al. (TPWRS 2022); Taylor/Lyapunov residual-error bounds (Scherpen 1993, Al'brekht 1961, Liu 2026); DAE reachability (Althoff-Krogh 2014).
Possibly new (unverified): the specific combination — an interval-certified IEEE-39 delayed-PLL counterexample where single-pole-preserving retunes do not compose, an exact rational interaction identity, an interaction-budget repair validated on nonlinear delayed events; the rank<=30 exact action space for SG->GFL replacement + PLL tuning; the exact finite-endpoint inertia-floor formula with graph secant resistance; the retraction/negative audit of a 94.2% candidate.
Weaknesses a reviewer will find: one bus pair; no comparison with published optimisers; the quadratic predictor also detects the counterexample; no maximum replacement or upper bound; floor not validated on IEEE-39; two loosely connected storylines.

## Open questions worth asking
1. Which single result is the most defensible and most novel for a student poster award?
2. What cheap, decisive experiment would make the pair-interaction finding systematic (all 45 bus pairs, preregistered; compare against a generic optimiser by number of full-model evaluations)?
3. Can the frequency-floor theorem be connected to the PLL/delay line (e.g. via effective damping d0 and primary response under PLL retuning)?
4. Is there 2024-2026 literature that already contains the pair-interaction determinant or the secant-resistance inertia floor?
