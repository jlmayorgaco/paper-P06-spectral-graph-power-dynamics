# PD39 PowerDynamics model and equilibrium audit

Status: qualified implementation scaffold; no PD39 stability or co-design campaign results are claimed here.

## Reference implementation

The PD39 baseline is the maintained IEEE-39 tutorial shipped with PowerDynamics.jl 5.0.0. The tutorial supplies 39 buses, 46 branches, 10 machine records, 19 load records, and the stock Sauer-Pai / AVR Type I / TGOV1 dynamic composition. Its CSV inputs are used verbatim through the package example; the existing Python IEEE-39 model in this repository is not substituted silently and none of its numerical results are reused as PD39 evidence.

The exact GFL replacement is the maintained library component
`PowerDynamics.Library.ComposableInverter.SimpleGFLDC`. It is compiled as a regular network bus, not as a custom model and not through an unrecorded parameter search. The frozen nominal settings are:

| Quantity | Value |
|---|---:|
| system base | 100 MVA |
| frequency base | 60 Hz |
| filter `Rf`, `Xf` | 0.01, 0.03 pu |
| PLL bandwidth / LPF | 5 Hz / 300 Hz |
| current-controller bandwidth | 600 Hz |
| DC voltage / capacitance | 2.5 pu / 1.25 pu·s² |
| DC voltage-controller bandwidth | 5 Hz |

PI gains are generated from those bandwidth values using the formulas in `src/pd39/model.jl`. They are frozen before the campaign. The stock model contains the PLL, L filter, current controller, DC-link state, and DC voltage PI state.

## Replacement semantics

For a candidate bus `i`, PD39 copies all compiled vertices and edges, obtains the original bus power-flow model, and constructs the GFL bus with `compile_bus(template; pf=original_pfmodel, vidx=i)`. Thus the replacement changes the dynamic device while retaining the official PF dispatch and voltage constraint at that bus. The slack bus is excluded from replacement.

The preregistered candidate set is the controlled-machine, non-slack set from the official bus table:

```text
{30, 32, 33, 34, 35, 36, 37, 38}
```

The official table identifies bus 31 as Slack. Bus 39 is an uncontrolled-machine-plus-load bus, not the reference bus.

## Qualification run

The reproducible qualification command is:

```text
julia --project=. scripts/pd39/00_model_qualification.jl
```

Recorded in `results/pd39/model_qualification/qualification.toml`:

| Check | Result |
|---|---:|
| Julia | 1.11.9 |
| buses / branches | 39 / 46 |
| baseline state dimension | 192 |
| single GFL replacement state dimension | 189 |
| baseline PF finite | true |
| baseline dynamic initialization finite | true |
| GFL replacement PF finite | true |
| GFL replacement dynamic initialization finite | true |

The bus-39 baseline voltage checkpoint is `u_r=1.014186196767751`, `u_i=-0.1797951594139631`. The bus-32 GFL qualification checkpoint is `u_r=0.982111980003339`. These are implementation checkpoints only, not scientific effect claims.

NetworkDynamics reports a deprecated `remove_conditions` keyword warning while constructing the Jacobian prototype. The warning is from the upstream API path; it does not invalidate the successful initialization, but it is retained as a provenance note and should be removed from the local wrapper when the upstream API is updated.

## Qualification boundary

This audit establishes that the selected maintained model can be assembled, solved, and initialized in the isolated PD39 worktree. It does not establish stability, a weak bus, a weak link, a reversal, or a feasible intervention design. Those quantities are governed by the preregistration in `docs/PD39_ROBUST_TRANSITION_PREREG.md` and must be generated only by the post-preregistration campaign scripts.

## Upstream references

- [PowerDynamics.jl repository](https://github.com/JuliaEnergy/PowerDynamics.jl)
- [IEEE-39 Part I: model creation](https://juliaenergy.github.io/PowerDynamics.jl/stable/generated/ieee39_part1/)
- [IEEE-39 Part II: initialization](https://juliaenergy.github.io/PowerDynamics.jl/latest/generated/ieee39_part2/)
- [Modeling concepts and custom component interface](https://juliaenergy.github.io/PowerDynamics.jl/stable/ModelingConcepts/)
