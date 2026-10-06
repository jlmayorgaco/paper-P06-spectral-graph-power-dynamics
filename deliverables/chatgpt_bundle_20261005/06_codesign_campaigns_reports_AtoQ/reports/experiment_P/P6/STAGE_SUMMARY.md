# P6 — Globality ledger

Status: **GLOBAL_GAP_OPEN**, for the nominal complete-spectrum objective only. All 1,024 support sets are enumerated in the ledger, but only one support has a joint continuous KKT candidate.

- Nominal spectral problem: valid global lower bound `L=0.0 MW`, from `P_i>0` and `ε_i≥0`; nominal spectral-feasible, independently PD-validated upper bound `U=1.124344477653483 MW`.
- Nominal spectral gap: `U−L=1.124344477653483 MW`; gap beyond a `0.01 MW` certificate tolerance: `1.114344477653483 MW`.
- The frozen candidate is **not** an upper bound for the combined robustness/transient-constrained problem: the declared 100 MW sustained step at bus 16 fails (`|Δf|=39.33340125578227 Hz`, RoCoF `1.072208443371192 Hz/s`, no settling within 60 s). P4 has one conditional point with a provisional Float64 interval lower bound above beta_req, but without outward rounding this is not a formal validated-arithmetic certificate; no robust optimum, joint KKT, or PD validation of that point exists. No combined-problem feasible incumbent or gap is available.
- The single-bus ExpP2 roots are conditional fixed-gain algebraic roots, not lower bounds for mixed multi-bus supports. The other 1,023 configurations remain open; the incumbent support also has no disconnected-branch completeness proof.
- GSP comparison: common gains (2 parameters), the ordered physical/noncommutative feature prefixes, the rank-5 full basis (10 coefficients) and the 20-gain reference all reproduce the same frozen uniform gain vector exactly; replacement loss is 0 MW. All PLLs remain local. This is a representation comparison at the incumbent, not evidence that the reduced family attains a different optimum.
- Generic-solver falsification was not run. The nominal global gap remains open and combined robust/transient optimization is incomplete; no global or robust-optimal claim follows.
- Candidate kept unchanged, hash `dbe727713acccec9b0802619f4b2e10763a0a88b9799ba24a74411ac88ff7046`. Time `26.205 s`.
