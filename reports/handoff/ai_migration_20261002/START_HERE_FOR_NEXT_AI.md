# Continuity packet for the next AI

Read this file first, then read HANDOFF_MASTER.md from top to bottom. The master dossier is the authoritative synthesis as of 2026-10-02; it distinguishes proven observations from hypotheses, records corrections, and links to the source artifacts. This package is an internal research handoff, not a camera-ready IEEE paper.

## Mission

Continue evaluating whether a defensible, network-aware method can compute a finite SG-to-GFL replacement boundary with nodewise PLL co-design, under explicitly declared DAE stability, frequency, RoCoF, voltage, current, DC-energy, actuator, disturbance, uncertainty, and recovery constraints. Keep the theory agnostic to particular buses and events; use IEEE-39 events only as finite validation scenarios.

## Read before changing or running anything

1. Check repository status and preserve the existing dirty tree. Do not clean, reset, checkout, bulk-format, or overwrite prior outputs.
2. Read HANDOFF_MASTER.md and ARTIFACT_MAP.json.
3. For any numerical claim, inspect the exact CSV/TOML/Julia files named in the map. Treat this handoff as a map, not a substitute for raw results.
4. Keep distinct the 90.785446% fixed-mix counterfactual, 90.047042% trajectory-gradient run, and 91.885% ExpQ2B best-found point. They use different candidates and contracts.
5. Review each campaign's terminal status before reusing it. FAIL_OPTIMIZATION, FAIL_TDS, failed robustness, or a failed gate must remain visible.

## Current evidence

The newest same-mix PowerDynamics experiment fixes 90.785446% GFL and tests six independent ±100 MW load steps at buses 8, 16, and 29. Co-designed PLL gains lower all six peak frequency deviations and change the exploratory 0.5-Hz frequency failures from 2/6 to 0/6; they also raise RoCoF in all six events, with maximum increasing 37.7%. Those six cases pass the exploratory voltage range. Current/DC hardware compliance, universal robustness, increased maximum replacement, and global optimality are not established.

The local linear recourse envelope allows only 0.001696 percentage points in its uniform replacement direction. This is a tangent diagnostic, not a nonlinear capacity bound. Analytic modal sensitivity agrees with centered finite differences to 2.61e-6 relative error. No useful second-order remainder or finite exclusion certificate exists yet.

## Research sequence

1. Reproduce the fixed-mix experiment without changing its frozen inputs; record Julia and package versions and hashes.
2. Establish a single explicit operating contract: dispatch, base values, devices, replacement vector, PLL gains, estimator/window, disturbance and uncertainty sets, hard limits, and simulation horizon.
3. Add diverse operating points and disturbances, plus declared current, DC-energy, and actuator limits. Re-test frequency/RoCoF tradeoffs.
4. Determine whether exact smooth active-branch DAE sensitivities and second derivatives can be bounded on a compact admissible box. Include algebraic elimination, mode switching, and estimator definitions.
5. Build a network-aware Hessian or interval enclosure for weighted active constraints. Compare its tightness and computational cost with dense, nodal, and generic interval alternatives.
6. Only then test a finite infeasibility/exclusion certificate for nearby replacement increases. Report certificate coverage, unresolved boxes, and a feasible candidate lower bound separately.
7. Audit novelty against current primary literature before describing the method as new. Standard Taylor bounds, KKT, sensitivity, graph Fourier transforms, and PLL co-design are not novelty by themselves.
8. Treat explicit digital latency as a separate extension: pure delay exp(-sT) differs from a first-order filter 1/(1+s tau). The current counterfactual includes neither a modeled PLL delay nor a sampled controller.

## Claim discipline

Use “in the six tested IEEE-39 events” for the observed performance result. Do not call the result globally optimal, robust to all disturbances, a maximum replacement boundary, or a validated Beyond Nodal Damping theorem. Distinguish a candidate feasibility lower bound from a global capacity upper bound. Do not infer device-limit compliance where no declared hard limit exists.

## Output expected from the next researcher

Maintain a dated decision log. For every campaign, save its code, immutable input contract, environment/version record, machine-readable outputs, figures, final status, and a short claim audit. Do not replace old results; create a new dated campaign directory. Update the master dossier and artifact map only after checking raw outputs.

Repository snapshot recorded in the master dossier: branch research/expQ2B-secure-optimum at HEAD d0fecb3264aeb855ab6700fc6ef66120d4b6c33c, with extensive working-tree modifications. The current live Git state may be newer; inspect it rather than assuming this snapshot is still current.
