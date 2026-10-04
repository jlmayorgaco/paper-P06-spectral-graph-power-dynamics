# Reproduction

Run from the repository root with Python3.13.9 and Julia1.11.9. Python versions
are recorded in BASELINE_MANIFEST.json; Julia's existing Project/Manifest and
package versions are preserved. The bundle contains required repository source,
exported matrices, design inputs and historical evidence used here. Install
the declared packages into your own environment; the archive is not an OS image.

```powershell
python -m pip install -r experiments/network_moment_codesign_20261004/requirements.txt
julia --project=experiments/physical_collective_damping_20261003/source_snapshot -e 'using Pkg; Pkg.instantiate()'
python experiments/network_moment_codesign_20261004/moments.py
python experiments/network_moment_codesign_20261004/verify_moment_balls.py
```

The input export can be independently regenerated with
`julia --project=experiments/physical_collective_damping_20261003/source_snapshot experiments/network_moment_codesign_20261004/export_inputs.jl`.
It derives the unit load-parameter inputs and complete39-bus phase outputs by
automatic differentiation and a78-dimensional algebraic solve. Input hashes
must match before accepting downstream results. Never replace a frozen input
with a materially changed model to make a check pass.

To reproduce the final numerical controller campaign:

```powershell
python experiments/network_moment_codesign_20261004/codesign_globalized.py
python experiments/network_moment_codesign_20261004/provenance.py
julia --project=experiments/physical_collective_damping_20261003/source_snapshot experiments/network_moment_codesign_20261004/validate_julia.jl
```

These commands write the experiment's result files. For a clean audit, unpack
the complete ZIP into a fresh folder and keep the distributed reference bundle
unchanged. `codesign.py`, `codesign_multimode.py` and their version snapshots
document failed/intermediate methods; they are not the final controller solver.
The code's input locks reject changed frozen definitions. Numerical root and
optimizer last digits can differ across BLAS and platform; compare tolerances,
feasibility and iteration status, not byte-identical trajectories.

`finalize.py` regenerates reports, figures and the section in the existing
THEORY.tex using THEORY_before_moment.tex as its fixed pre-edit source. It does
not recompute experiments. Compilation reads that same source from repository
root so its existing figure paths resolve. No replacement document is opened.

The DDE spectral oracle uses the full exponential, gauge deflation and a
computed bounding region. Its complete-region count is numerical, not Arb
certified. Nonlinear validation uses the existing method-of-steps/adaptive
Rosenbrock solver with uniform41ms delay and frozen five-event metrics.
Heterogeneous cases have catalog checks only; no heterogeneous nonlinear
simulation is implied. The moment ball check certifies an explicitly
symmetry-restored binary export, not uncertain hardware or plant parameters.
