# PD39 255+1 mechanism-validation execution plan

## Frozen starting point

1. Verify `902403cf` and work only in the clean worktree
   `C:\w\pd39v`.
2. Keep branch `research/pd39-255plus1-mechanism-validation` active.
3. Read frozen discovery/confirmatory artifacts without modifying them.
4. Commit the preregistration, claim matrix, and this plan before any new
   numerical process.

## Execution order

### A. Numerical truth audit

* portfolios: V8 plus `missing_30`, `missing_32`, `missing_33`, `missing_34`,
  `missing_35`, `missing_36`, `missing_37`, `missing_38`;
* scenarios: all nine rows in the frozen discovery scenario definition;
* tolerances: `1e-8`, `1e-10`, `1e-12` for the PF/initialization sweep;
* independent Jacobian: central finite differences of the network RHS;
* eigensolvers: reduced `eigen`, `eigvals`, and generalized `(A,M)` when the
  descriptor matrices are valid;
* output: one auditable row per portfolio/scenario/tolerance/method.

### B. High-PLL nonlinear TDS

* four `pll=1.2` corners;
* V8, best 7/8 by frozen `m_9`, worst 7/8 by frozen `m_9`;
* 12 total simulations, +1% P/Q pulse, 100 ms, 1–20 s;
* output trace and summary for every run.

### C. Complete holdout census

* input: existing 24-condition file, unchanged;
* portfolios: every 256 bit-mask portfolio;
* attempts: exactly 6144 unique `(portfolio, condition)` keys;
* output: full row-level census plus condition-wise `H0` and `H0.05` tables.

### D. Mechanism

* aggregate-matched 6/8 and 7/8 pairs selected from static metadata only;
* all eight 7/8 predecessors and V8 tracked at nominal and high-PLL cases;
* common bus-voltage mode shape, MAC, participation, alpha/frequency;
* homotopy only if the installed model exposes a documented continuous SG↔GFL
  interpolation; otherwise write `BLOCKED`.

### E. Exact closure

* inspect source/API for exact definitions of `K`, `D`, and `Q`;
* construct only if dimensions, units, and semantics are explicit;
* test determinant identity and `Q→-1` boundary only after feasibility;
* otherwise report `BLOCKED` with no proxy.

## Excluded work

Do not run weak-node/link analyses, 128-direction structured radii, repair
optimization, planners, co-design, IEEE-68, EMT, or new controller search.

## Final artifacts

The campaign must create the final PDF report, ChatGPT handoff, headline JSON,
master claims CSV, and upload ZIP.  The final response must explicitly include
the exact line `FINAL CASE: A / B / C` and the requested pass/fail/block status
diff.

