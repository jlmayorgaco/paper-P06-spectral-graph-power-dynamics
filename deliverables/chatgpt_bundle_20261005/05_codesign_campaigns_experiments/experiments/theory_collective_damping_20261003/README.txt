THEORY FIRST — 2026-10-03

Open THEORY.tex in the built-in LaTeX editor for the complete equations,
assumptions, proofs, counterexample, literature boundary and claim ledger.

Primary candidate: exact single-channel latency redistribution theorem
at physical SG torque-speed ports. Noncollinear actuation/measurement
paths imply an indefinite rank-two change of harmonic damping.

Complement: a graph Green-function correction to the aggregate frequency
moment, and a signed-area theorem forcing overshoot in an explicit
reduced model. This reduced model is NOT yet justified for the full
IEEE-39 lossy PLL/DC model.

Actual PLL correction: closed phase-tracking delay enters the Taylor
series first at cubic order, subject to regularity of later reductions.

Run algebra checks:
    python experiments/theory_collective_damping_20261003/verify_symbolic.py

Requires SymPy. No Julia run, nonlinear simulation, optimizer, delay
sweep or eigenvalue continuation is performed.

SYMBOLIC_CHECKS.json records exact algebra checks, not empirical IEEE-39
results. PROVENANCE.json records inspected source hashes and dirty state.

No new rho*, Kp*, Ki*, maximum replacement or global optimality gap
is claimed. Bibliographic originality of the candidate theorem remains
unestablished. All prior frozen files are preserved.

THEORY-ONLY RETUNING EXTENSION
The same open THEORY.tex now also contains:
  - exact one-site pole assignment through the full-grid PLL return relation;
  - real-gain feasibility as a complex parallelogram;
  - closed-form gain transport between prescribed heterogeneous latencies;
  - the exact modal spillover identity and its single-site limitation;
  - a gain-bounded maximum replacement step for first-order constraints;
  - a conditional complete-DDE-spectrum preservation test.

Run its additional algebra checks:
    python experiments/theory_collective_damping_20261003/verify_retuning_symbolic.py

RETUNING_SYMBOLIC_CHECKS.json: thirteen exact algebra checks, no grid experiments.
RETUNING_PROVENANCE.json records the revised source hash. PROVENANCE.json
and SYMBOLIC_CHECKS.json are preserved as records of the initial version;
their THEORY.tex hash intentionally refers to that earlier source.

Pole assignment and convex separation are established methods. The proposed
research contribution concerns the collective spillover mechanism and its
use in replacement decisions; bibliographic novelty remains open.

ANALYTICAL DESIGN CLOSURE
Section "Closing the design problem before simulation" in the same THEORY.tex
now supplies the all-PLL return, explicit gains for a collective eigenpair,
multi-pattern gain compatibility, and the minimum added-neighbor gain law.
It also states the limits of polynomial graph filters with heterogeneous
measurement delays and gives the gain-map differential for a replacement
predictor/corrector. No new controller or secure replacement is validated.

Run only its algebra verification:
    python experiments/theory_collective_damping_20261003/verify_closure_symbolic.py

CLOSURE_REVIEW.txt records the research contract, adversarial scope audit,
verified literature boundary and deferred validation gates.
CLOSURE_SYMBOLIC_CHECKS.json and CLOSURE_PROVENANCE.json record this revision.
Earlier provenance/check files remain unchanged and refer to their own versions.
