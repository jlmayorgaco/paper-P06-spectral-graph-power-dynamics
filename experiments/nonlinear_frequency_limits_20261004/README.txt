NONLINEAR FREQUENCY LIMITS — DEVELOPMENT PACKAGE

Main source remains the user's existing open document:
  experiments/theory_collective_damping_20261003/THEORY.tex
The reviewed section is copied into that file in place by
  integrate_manuscript.py
The pre-edit source is preserved here as THEORY_before_extension.tex.
This backup is not a replacement document or a new editor target.

THEORY
  FREQUENCY_LIMITS_SECTION.tex: full hypotheses, proofs, graph cuts,
    measurement window, memory, energy and retention-relaxation scope.
  MATH_REVIEW.txt: independent adversarial mathematical review.
  NOVELTY_AND_POSTER_REVIEW.txt: primary-source comparison and limitations.
  POSTER_ARGUMENT.txt: central claim, visual story and remaining gate.

EXECUTED CHECKS — NO NEW DYNAMIC SIMULATION
  python experiments/nonlinear_frequency_limits_20261004/verify_closed_form.py
Evaluates the exact feasible construction and checks the nonlinear equations,
integrals and two-node energy barrier. Also solves24 static endpoint cases
on a six-node illustration to check the maximum principle and graph formulas.
Outputs: CLOSED_FORM_RESULTS.json, STATIC_GRAPH_CHECKS.csv, PNG/SVG figure.
These are illustrative mathematics, not new IEEE39 trajectories.

ARCHIVED IEEE39 AUDIT
  python experiments/nonlinear_frequency_limits_20261004/audit_archived_ieee39_balance.py
Reads seven frozen CSVs and their source/parameter contracts. Reconstructs
kinetic, filter and DC quantities, losses, torque and governor balances.
IEEE39_BALANCE_REVIEW.txt explains the sign convention and physical gate.
Its JSON contains input hashes, parameters and versions; CSV files contain
derived quantities, not re-simulations. The old source/results were not edited.

MANUSCRIPT
  python experiments/nonlinear_frequency_limits_20261004/integrate_manuscript.py
Only the existing THEORY.tex is updated. The section, introduction pointer
and verified references are integrated; older numerical records remain intact.
Native compile was attempted before and after editing. It fails at host
initialization ('Unable to find standard directories for platform'), before
TeX source diagnostics. This package does not claim a compiled updated PDF.

STATUS
EXACT_IDENTITY: fixed-output endpoint areas, memory/filter corrections.
PROVED_REDUCED_MODEL: all-bus necessary floor, graph/cut forms, sharp2node
infimum and explicit feasible control without actuator bounds.
NUMERICALLY_VALIDATED: static/algebra checks and archived model balance.
BLOCKED: physical IEEE39 energy transfer, global SG/GFL optimality,
publication priority, updated PDF compilation and final engineering poster.

No commits, pushes, installations or changes to frozen model/results.
