# Model provenance

Git HEAD: `d0fecb3264aeb855ab6700fc6ef66120d4b6c33c`. Repository already dirty; full status and
exact source hashes are in EXPERIMENT_MANIFEST.json. This experiment writes only
inside its own directory. Prior frozen designs/results are copied, not modified.

The model uses the repository's exact algebraic network elimination and full
SG/AVR/governor and GFL/DC/current/PLL differential equations. Physical-supply
DC convention. Pure delay is inserted only at the PLL voltage-error detector.
The original PLL frequency low-pass state remains, independently of pure delay.
No physical line dynamics or current limiter absent from the source model is
added or claimed. The ten retained SG ratings define physical torque-speed ports.

Python: 3.13.9; NumPy 2.3.5; SciPy 1.16.3.
Julia manifest version: 1.11.9. Package versions: {'CSV': '0.10.17', 'DataFrames': '1.8.2', 'ForwardDiff': '1.4.6', 'NetworkDynamics': '1.3.0', 'OrdinaryDiffEqRosenbrock': '2.7.1', 'PowerDynamics': '5.0.0', 'SciMLBase': '3.53.2'}.

`source_snapshot/` includes the exact minimal repository dependencies and inputs.
The downloadable zip lays those out at the repository-relative paths required
by the code. Project/Manifest files and original source paths are retained.
