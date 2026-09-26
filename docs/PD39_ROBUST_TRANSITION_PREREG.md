# Preregistration: robust SG→IBR transition co-design on PowerDynamics IEEE-39

Status: preregistration for the PD39 campaign. This document is frozen before campaign-level stability, weakness, or design results. Any deviation must be appended with date, reason, affected claims, and whether the primary result changes.

## Scope and separation from prior work

PD39 is an independent benchmark line. Prior IEEE-39 and IEEE-68 Python/custom-model numbers in this repository may motivate questions or define comparison language, but they are not PD39 measurements. The TX4 frozen model, its negative results, and its author decisions remain untouched. No claim from PD39 may be backfilled into the TX4 evidence set.

## Frozen model and candidate set

- Baseline: the maintained PowerDynamics.jl IEEE-39 tutorial CSV model, with 39 buses and 46 branches.
- Dynamic SGs: the stock Sauer-Pai / AVR Type I / TGOV1 composition from that tutorial.
- Exact GFL: `PowerDynamics.Library.ComposableInverter.SimpleGFLDC`.
- Nominal GFL parameters: the values in `docs/PD39_POWERDYNAMICS_MODEL_AUDIT.md` and `src/pd39/model.jl`.
- Replacement: one controlled-machine bus is reconstructed as a GFL bus while retaining its original PF model. No dispatch redispatch, slack-bus replacement, topology repair, controller retuning, or candidate-specific parameter fitting is allowed in the primary campaign.
- Slack bus: 31. It is excluded from intervention.
- Candidate SG buses: `{30, 32, 33, 34, 35, 36, 37, 38}`.
- Portfolio set: all `2^8=256` subsets, enumerated by increasing integer bit mask; no random subset sampling.

## Equilibrium and stability definitions

For every network/scenario, the script first solves the PowerDynamics power flow and then calls `initialize_from_pf!`. An equilibrium is usable only if the PF state and dynamic state are finite and the fixed-point check passes at tolerance `1e-8`.

For a usable equilibrium, `jacobian_eigenvals` supplies the reduced DAE spectrum. Eigenvalues with `|λ| <= 1e-8` are recorded as numerical/gauge modes and excluded from the engineering margin. The raw spectrum and the excluded modes remain in the output. With the remaining eigenvalues,

```text
m_dyn = - max Re(λ_nontrivial)
```

The nominal linear-stability predicate is `m_dyn > 1e-8 s^-1`. The robust design threshold is fixed at `m_target = 0.05 s^-1`; this is a design requirement, not an observed baseline property.

## Robustness set

The primary robust screen contains nine controller scenarios: the nominal case and all eight corners of independent ±20% factors on PLL bandwidth, filter reactance, and current-controller bandwidth. PLL `Kp` and `Ki` are regenerated consistently from the perturbed bandwidth; the other parameters remain frozen. The uncertainty grid is:

```text
{1.0} ∪ {0.8, 1.2}PLL × {0.8, 1.2}Xf × {0.8, 1.2}CC
```

A portfolio is robustly feasible only if every one of its nine scenarios has a finite usable equilibrium, `m_dyn >= 0.05 s^-1`, and the non-gauge spectrum is linearly stable.

## Weak nodes and links

Static rankings are computed before dynamic results:

- node proxy: inverse incident electrical strength, `1 / Σ_e 1/|X_e|`;
- link proxy: absolute series reactance `|X_e|`.

The dynamic node score is single-replacement margin loss, `m_baseline - m_i`, and is computed only from the preregistered candidate buses. Dynamic link analysis is a separate diagnostic using one-at-a-time removal of each of the 46 official branches, with PF re-solved and no redispatch. It preserves the same equilibrium/stability checks and must not be selected after inspecting the node result. If a branch-removal case cannot be qualified, it is reported as not estimable rather than replaced by a post-hoc proxy.

The primary topology/dynamics comparison is rank association (Spearman ρ) between the static score and the dynamic score, with ties handled by the implementation's fixed ranking routine. A low or sign-reversed association is evidence about contextual dynamic weakness, not proof of a universal weakness index.

## Minimum-intervention design

The design objective is lexicographic:

1. minimize the number of SG→GFL replacements;
2. among robustly feasible portfolios at that cardinality, maximize the worst-case `m_dyn`;
3. use portfolio bit-mask order as the deterministic final tie-break.

There is no feasibility fallback. If no portfolio satisfies the robust threshold, the primary design result is “no feasible design in the preregistered candidate set and uncertainty box.” Any relaxed threshold becomes a new exploratory analysis and cannot replace the primary result.

## Holdout and stopping rules

Linear spectra are the primary screen. Time-domain simulations are holdout checks only after the campaign table is frozen; their disturbance, duration, solver, and acceptance criteria must be recorded before running them. A simulation failure is not silently converted to a stable label.

The campaign stops after all 256 portfolios × 9 scenarios are attempted, including failed equilibria and unstable cases. No portfolio or scenario may be omitted because it is inconvenient, slow, or unfavorable. Failed implementation runs are retained with an error status and are separated from physical instability.

## Planned outputs

- raw portfolio/scenario table with equilibrium, spectrum, margin, and error provenance;
- static node and link tables;
- single-replacement dynamic node table;
- robust-feasibility/design summary;
- deviation log and claims matrix;
- no claim of universal weakness, irreducible interaction, or causal mechanism without an independently qualified test.

## Reproducibility

The Julia environment is pinned by `Project.toml` and `Manifest.toml` in the PD39 worktree. The model audit must pass before `scripts/pd39/01_baseline.jl` or any portfolio campaign script is run. Results are written under `results/pd39/` and must include Julia/package provenance and the exact scenario identifiers.
