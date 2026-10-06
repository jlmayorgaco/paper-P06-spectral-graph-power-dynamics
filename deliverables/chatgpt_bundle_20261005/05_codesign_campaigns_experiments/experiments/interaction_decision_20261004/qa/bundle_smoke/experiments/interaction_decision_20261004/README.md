# Finite collective interaction and a PLL compensation rule

Read REPORT_ES.md and POSTER_CLAIMS.md first. The new theory is incorporated
in place into experiments/theory_collective_damping_20261003/THEORY.tex,
Section9. The38-page compiled PDF in this directory includes the preserved
earlier derivations. This is a research manuscript, not a submission-ready claim.

The executed result is a certified finite modal decision reversal and a
corrected full-model numerically feasible design. It is NOT a maximum of GFL
penetration. Fixed delays were not optimized. Existing results were preserved.

## Reproduction

First run `python experiments/interaction_decision_20261004/restore_dependencies.py`
to restore missing archived source/data paths from DEPENDENCIES.zip. It refuses
to overwrite differing existing files. Install requirements.txt in your chosen
Python environment. On a fresh Julia installation, instantiate the archived
Project/Manifest environment before running the Julia command below.

Run from the repository root in a fresh copy to avoid overwriting recorded
campaign outputs. Python3.13.9, numpy2.3.5, scipy1.16.3, pandas2.3.3,
matplotlib3.10.6, python-flint0.8.0. The original host finds flint under
experiments/regional_paper_closure_20261003/vendor; normal installed flint also
works when that path is absent. No Padé approximation is used.

1. python experiments/interaction_decision_20261004/discover.py
2. python experiments/interaction_decision_20261004/compensated.py
3. python experiments/interaction_decision_20261004/damping_only.py
4. python experiments/interaction_decision_20261004/interaction_contour.py
5. python experiments/interaction_decision_20261004/certify_roots.py
6. python experiments/interaction_decision_20261004/verify_common.py
7. python experiments/interaction_decision_20261004/close_design.py
8. python experiments/interaction_decision_20261004/compare_predictors.py
9. python experiments/interaction_decision_20261004/policy_path.py

certify_roots.py certifies the five original endpoint designs. The additional
corrected endpoint is certified by calling certify('corrected',p,root) using
BUDGET_RESULT.json; policy_path.py certifies its two endpoints. See finalize.py
for the end-to-end evidence audit and artifact manifest generation. The
CODE_AND_RESULTS complete ZIP retains directory paths relative to the repository
root. Its actual filename is INTERACTION_DECISION_20261004_COMPLETE.zip.

Independent full model and nonlinear events (~10min on the original host):

    julia --project=experiments/physical_collective_damping_20261003/source_snapshot experiments/interaction_decision_20261004/validate_julia.jl

Build tables, scientific plots and the existing manuscript from saved evidence:

    python experiments/interaction_decision_20261004/build_report.py

The native document compiler was attempted and returned 'Unable to find
standard directories for platform'. The already installed MiKTeX compiled the
same source successfully. Run the following twice from the repository root:

    pdflatex -interaction=nonstopmode -halt-on-error -output-directory=experiments/interaction_decision_20261004 -jobname=Beyond_Nodal_Damping_Interaction_Decision_20261004 experiments/theory_collective_damping_20261003/THEORY.tex

All initial discovery cases, adaptive addenda, failures, certificate panels,
gain vectors, full spectral counts and five event trajectories are retained.
This campaign was adaptive research, not an independently held-out prediction
benchmark. READ REVIEW_RECORD.md for deviations and unresolved solver cases.
