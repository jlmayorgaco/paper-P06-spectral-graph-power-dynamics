# Experiment Q — IEEE-39 exact action-space SG→GFL co-design

This is a new, isolated experiment. It reads the frozen repository model and previous experiment outputs as inputs, and writes every new artifact here. It does not modify prior results.

The first gate is a fully feasible zero-delay point under the five frozen physical events: bus 8 at −100 MW; bus 16 at ±100 MW; and bus 29 at ±100 MW. A point is not feasible unless the complete finite spectrum, all five full PowerDynamics simulations, frequency and RoCoF limits, voltage guard, gain bounds, and the declared SG actuator margin pass.

Run the zero-delay validation from the repository root:

```powershell
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q0_validate_five_events.jl experiments/ieee39_exact_action_space_codesign_20261002/seed_uniform_875.toml
```

`INCOMPLETE` is not a pass. This directory's `STATUS.md` is the current gate ledger. No global-optimality claim is made.
