THEORY FIRST — REGIONAL REPLACEMENT BOUND

Read REPORT_ES.txt, STATUS.json, and the new section in the existing THEORY.tex.
SECTION.tex is a source fragment integrated in place, not a replacement document.

Proof:
  root-cluster log-determinant identity;
  finite contour remainder;
  exact affine descriptor lift for the inspected sharing contract;
  40-channel repeated-diagonal parameter perturbation;
  positive-matrix box bound and continuous-panel error expression;
  regional optimistic LP, gain-eliminated dual and exclusion witness;
  trace closed-walk interpretation with explicit limitations.

Executed:
  checks.py -- exact algebra, two-state DDE fixture, 3 IEEE-39 matrix identities.
  check_descriptor.py -- same points/box, affine descriptor and sharper bounds.
  integrate.py -- generated report and in-place document update.

Run from the repository root with the existing Python environment:
  python -B experiments/theory_collective_damping_20261003/regional_certificate_20261003/checks.py
  python -B experiments/theory_collective_damping_20261003/regional_certificate_20261003/check_descriptor.py
  python -B experiments/theory_collective_damping_20261003/regional_certificate_20261003/integrate.py

Run numerical reproduction in a copy to preserve the executed outputs.
integrate.py refuses to overwrite manual document edits made after integration.
The native editor compiler must be used after updating the open source.

The protocol addendum documents the analytical refinement after the naive
scalar envelope failed. The original pointwise results remain in TABLE_02.
Neither TABLE_02 nor TABLE_03 verifies an entire contour.

No power-system trajectories, optimization, new feasible design or IEEE-39
regional maximum have been computed here. No result is claimed globally optimal.
