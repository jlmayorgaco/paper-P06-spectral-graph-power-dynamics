# Reproduce Experiment G

Run from the repository root on the pinned Julia 1.11 environment:

    julia --project=. --startup-file=no experiments/bnd_expG/run_experiment_G.jl

After Z_G_FINAL.toml exists, reproduce the analytical tables while verifying and preserving that frozen candidate:

    julia --project=. --startup-file=no experiments/bnd_expG/run_experiment_G.jl --resume-frozen
    python experiments/bnd_expG/render_figures_G.py
    julia --project=. --startup-file=no test/bnd_expG/runtests.jl
    julia --project=. --startup-file=no experiments/bnd_expG/validate_powerdynamics_G.jl

The design runner reads frozen D/E input artifacts but writes only under experiments/bnd_expG and reports/experiment_G. It does not import PowerDynamics. The separate validator refuses to run unless Z_G_FINAL.toml and its SHA-256 sidecar verify.

The all-GFL endpoint has a defective gauge/physical-zero pair in the full state matrix. ExpG deflates the exact global-angle gauge and reports the remaining physical zero separately. The scalar authority LP is a seed; the final analytical candidate is checked against the complete finite spectrum and robust conditions.

The robust radius is dimensionless relative to the declared 1 s^-1 additive full-block state-matrix scale. It is not a physical percentage standard. No unequilibrated load-setpoint variation is asserted.
