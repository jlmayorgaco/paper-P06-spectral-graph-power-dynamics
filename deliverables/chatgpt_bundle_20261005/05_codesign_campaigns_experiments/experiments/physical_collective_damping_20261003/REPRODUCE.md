# Reproduction

Run from repository root (or root of the extracted bundle). Julia 1.11 with the
included Project.toml/Manifest.toml and Python with numpy, scipy, pandas,
matplotlib are required. See EXPERIMENT_MANIFEST.json for executed versions.

```powershell
julia --project=. experiments/physical_collective_damping_20261003/export_model.jl
python experiments/physical_collective_damping_20261003/analyze_ports.py
julia --project=. experiments/physical_collective_damping_20261003/nonlinear_steps.jl
python experiments/physical_collective_damping_20261003/action_attribution.py
python experiments/physical_collective_damping_20261003/finalize.py
```

`export_model.jl` checks matching hashes when a frozen design copy already exists
and refuses to replace a different design. No original input should ever be
modified. Other scripts replace only outputs
inside this new experiment directory; archive the bundle before a rerun.

For an extracted bundle without Git, run the computational steps and omit the
provenance/package step in finalize.py, or use `python finalize.py --no-git`
(see script support). Frozen provenance from the original run remains in the bundle.

The spectral calculation refines roots of the exact exponential characteristic,
using deterministic frequency seeds, frozen critical roots and adaptive
predictor/corrector continuation. It is not an exhaustive infinite-spectrum proof.
The nonlinear histories are small eigenmode perturbations, not load-fault cases.

Source dependencies are mirrored under source_snapshot/ in the original output
directory; the zip places them at the correct repository-relative locations.
