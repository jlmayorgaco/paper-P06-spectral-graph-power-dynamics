# Reproduce Experiment H

Run from the repository root on branch `research/expH-jordan-mixed-codesign`.

1. Run the analytical, PowerDynamics-free design and its preblind/final freezes:

   ```powershell
   julia --project experiments/bnd_expH/run_experiment_H.jl
   ```

   Repeated runs verify the frozen candidate hashes and preserve already-written post-freeze H20/H21 results.

2. Run independent post-freeze PowerDynamics and nonlinear TDS validation. This script verifies the exact final candidate SHA before loading PowerDynamics and does not change candidate values:

   ```powershell
   julia --project experiments/bnd_expH/validate_powerdynamics_H.jl
   ```

3. Render all report figures from the H00–H21 CSV tables:

   ```powershell
   python experiments/bnd_expH/render_figures_H.py
   ```

4. Check the frozen evidence and scientific gates:

   ```powershell
   julia --project -e 'using Test; include("test/bnd_expH/runtests.jl")'
   ```

Expected independent validation result is `PARTIAL`: the equilibrium and TDS pulse scaling pass, while PowerDynamics spectral abscissa is approximately \(-0.001326\,s^{-1}\), outside the requested \(-0.05\,s^{-1}\) margin and inconsistent with the analytical \(-0.064796\,s^{-1}\) value. The experiment must not be described as fully validated.
