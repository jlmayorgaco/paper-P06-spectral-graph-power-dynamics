BEYOND NODAL DAMPING — REGIONAL PAPER CLOSURE
3 October 2026

Start with:
  ../theory_collective_damping_20261003/THEORY.tex — the existing open manuscript
  STATUS.json — machine-readable evidence gates
  POSTER_CLAIMS.txt — the permitted claims
  REVIEW_RECORD.txt — three independent internal reviews and remaining gates
  LITERATURE_COMPARISON.txt — primary-source overlap and contribution boundary

NEW MATHEMATICS
1. Fixed-chart affine descriptor for replacement and local PLL gains.
2. Finite contour-cluster identity with free modal patterns.
3. Exact cluster Hessian resolving mixed replacement/gain actions.
4. Closed-walk remainder and a rigorous budget for dropping inter-site terms.
5. Necessary regional LP and quadratic semidefinite outer relaxation.

These are carefully scoped derivations, not claims that descriptor lifting,
contour methods or convex relaxation were invented here. Their application must
still yield an informative physical replacement decision to support the strongest
paper claim. Current 87.6% evidence is a constructed candidate; it is not a
demonstrated physical maximum. The small verification box is limited by its own
rho ceiling and cannot establish irreducible SG support.

NEW EXECUTED CHECKS
TABLE_01_WALK_REMAINDER.csv — three frozen IEEE-39 algebraic frequencies.
TABLE_02_QUADRATIC_TOY.csv — exact/contour quadratic derivatives and finite checks.
PANELS.json / CONTOUR_RESULT.json — continuous complex-ball contour enclosures.
No nonlinear power-grid trajectories or replacement optimization were run during
this closure. The previous nonlinear evidence is archived separately and retains
its original scope and hashes.

REPRODUCTION
Use a separate copy before rerunning commands: outputs are regenerated in place.
Python package versions are pinned in requirements.txt. Install them in an
isolated environment, or install python-flint task-locally as:
  python -m pip install --no-deps --target experiments/regional_paper_closure_20261003/vendor python-flint==0.8.0
The other packages must match requirements.txt.

From repository root:
  python -B experiments/regional_paper_closure_20261003/theory_checks.py
  python -B experiments/regional_paper_closure_20261003/verify_contour.py
  python -B experiments/regional_paper_closure_20261003/build_manuscript.py
  python -B experiments/regional_paper_closure_20261003/make_figures.py

The contour run uses complex balls at 128-bit precision and may take several
minutes. It verifies one declared contour, not the entire DDE spectrum.
The input CSV binary64 values are treated as the exact linear model. This does
not enclose equilibrium uncertainty or error relative to a physical grid.
Ordinary sampled poles or midpoint arithmetic do not substitute for its gates.

Open THEORY.tex in the native LaTeX editor to compile. The calling environment
currently reports "Unable to find standard directories for platform"; therefore
PDF compilation/visual verification is not confirmed. No separate PDF is built.

HISTORICAL NONLINEAR EVIDENCE
evidence/all_pll_gain_map_validation_20261003_frozen.zip contains the complete
preserved experiment outputs and scripts, including the five saved trajectories.
evidence/HISTORICAL_REPRO_SOURCE.zip contains the model equations, network data,
direct source dependencies and frozen Julia Project/Manifest files.
To reproduce that earlier experiment, extract both archives into a NEW working
directory and follow its experiments/all_pll_gain_map_validation_20261003/README.txt.
Use Julia 1.11.9 with the recorded environment. This closure does not reinstall
Julia dependencies or rerun that experiment. Source dependencies were audited
and hashed; a fresh-environment Julia rerun remains a separate verification.
evidence/regional_certificate_20261003_frozen.zip preserves the preceding
conditional theorem/algebra experiment.

PRESERVATION
inputs/ contains byte-identical direct input copies for portable new checks.
BUNDLED_INPUTS.json and HISTORICAL_SOURCE_MANIFEST.json record source hashes.
THEORY_before_closure.tex preserves the document before this turn.
PROTOCOL.txt was frozen before new checks; AUTHORIZATION_UPDATE.txt records the
subsequent user request to commit/push and the endpoint-arithmetic repair.
FINAL_MANIFEST.json records the closure artifacts. vendor/ is deliberately
excluded from Git; dependency versions, not platform binaries, are published.
