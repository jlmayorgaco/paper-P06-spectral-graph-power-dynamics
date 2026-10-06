# Figure Design Contract: H4 Spectral Power-Dynamics Mechanism

**Status:** design freeze before the next scientific run  
**Scope:** poster figures, supporting evidence strip, and the data required to regenerate them  
**Current gate:** M1A is complete but `BLOCKED_M1_STRICT`; M2 must not start until the matched-load-block discrepancy is repaired.

## 1. Poster-level composition

The poster uses four primary figures and one narrow evidence strip:

```text
+----------------------+----------------------+
| F1 Stability Cliff   | F2 Local vs collective|
| 16 portfolios        | closure / singular    |
| inclusion edges      | values                |
+----------------------+----------------------+
| F3 Feedback map      | F4 GFL/GFM rescue     |
| only after M1 passes | 4 x 4 technology grid |
+----------------------+----------------------+
| Evidence strip: Python--Julia parity | TDS, only when valid |
+------------------------------------------------+
```

The visual order is causal rather than chronological:

1. F1 establishes the stability cliff over the H4 inclusion lattice.
2. F2 separates physical local closure from collective return closure.
3. F3 tests the reduced feedback mechanism in the `(g, eta)` plane.
4. F4 tests whether the mechanism survives technology substitutions.
5. The evidence strip documents reproducibility and cross-code agreement.

No figure may be used to imply a result that its readiness gate has not cleared.

## 2. Figure specifications

### F1 — Stability Cliff: H4 inclusion lattice

**Question:** Is the H4 instability caused by a single critical inclusion, or does it require the complete selected portfolio?

**Panel design**

- Main panel: 16 portfolio points, x-axis = portfolio cardinality `|H|`, y-axis = dominant real eigenvalue `alpha(H)` in `s^-1`.
- Draw inclusion edges from each proper subset to each one-element extension. Use thin neutral edges; do not connect portfolios solely because they have equal cardinality.
- Use one semantic color for stable points (`alpha < 0`), one for the boundary (`|alpha| <= tolerance`), and a high-contrast accent for H4 (`alpha > 0`).
- Add a horizontal zero line and a compact H4 annotation with `alpha` and `f = Im(lambda)/(2 pi)`.
- If labels become crowded, label only H4, its immediate predecessors, and the most negative/positive point; retain the complete identifiers in the data table.

**Required data**

```text
portfolio_id, portfolio_mask, cardinality, parent_ids, child_ids,
alpha, dominant_lambda_real, dominant_lambda_imag, dominant_frequency_hz,
spectral_abscissa_tolerance, status, equilibrium_hash, model_hash,
solver_backend, run_id
```

**Acceptance criteria**

- Exactly 16 portfolios for the H4 lattice.
- Full spectra are retained for every portfolio, not only the dominant eigenvalue.
- The H4 point is reproduced by the same equilibrium and operator contract used by the baseline campaign.
- The plot caption states whether the lattice is Python-only or cross-code reconciled.

**Current readiness:** design ready; data should be regenerated or collected from the exact GFL11 campaign before poster export. Do not infer the full 16-case lattice from a selected-portfolio CSV.

### F2 — Local versus collective closure

**Question:** Does the local physical closure become singular before the collective return closure does?

**Panel design**

- x-axis: coupling/feedback parameter `g`.
- y-axis: minimum singular value on a logarithmic scale.
- Solid curve: physical local factor `sigma_min(I + M_ii)` for each monitored bus/port, with the limiting port highlighted.
- Dashed curve: collective factor `sigma_min(I + Q_H)` for the selected H4 block.
- Mark the observed boundary `g*` with a vertical line and annotate the corresponding `alpha = 0` crossing.
- If multiple local ports are shown, keep a fixed port ordering and use small multiples or a direct label; do not rely on an unexplained legend.

**Important normalization rule**

The local curve must use the physical pre-normalized operator `I + M_ii`. A diagonal block of the normalized return matrix `Q` is not a substitute for the physical local closure. The data source must identify which convention was used.

**Required data**

```text
g, port_id, local_sigma_min_physical, collective_sigma_min,
local_operator_hash, collective_operator_hash, boundary_g_star,
alpha_at_g, equilibrium_hash, model_hash, run_id
```

**Acceptance criteria**

- Physical and collective curves are calculated from the same equilibrium, port ordering, and load block.
- The figure includes the operator convention in its caption.
- The baseline physical-local factors are available in `TX4_PHYSICAL_LOCAL_FACTORS.csv`; the older normalized-diagonal sweep is not a valid replacement.

**Current readiness:** design and baseline evidence are available. Treat this as a mechanism figure, not as proof of the feedback reduction until F3 clears its M1 gate.

### F3 — Feedback stability map

**Question:** Does the reduced feedback model predict the full-network stability boundary over the `(g, eta)` plane?

**Gate:** `PENDING_M1_REPAIR`. This figure must not be generated or interpreted until the matched-load-block discrepancy in M1A is resolved and M1 strict passes.

**Panel design**

- Main panel: heatmap of `Re(lambda_feedback(g, eta))` over a fixed, versioned rectangular grid.
- Overlay the zero contour `Re(lambda_feedback) = 0` with a thick contour line.
- Mark the full-network boundary and the nominal H4 operating point.
- Add a small inset with the tracked complex root locus along the `eta` slice used for the headline claim.
- Keep the color scale symmetric around zero and report the numerical range in the metadata.

**Required data**

```text
g, eta, feedback_lambda_real, feedback_lambda_imag,
full_network_lambda_real, full_network_lambda_imag,
root_id, root_tracking_method, boundary_residual,
reduction_residual, m1_status, MAC_to_H4, eigenvalue_gap, mode_family_id,
equilibrium_hash, model_hash, run_id
```

**Acceptance criteria**

- M1 strict residual is below the declared tolerance at the base point and along the sweep.
- Root identity is tracked continuously; independently re-solved roots are retained for audit.
- The map includes the full-network comparison, not only the reduced model.
- The figure is marked `BLOCKED` in metadata if the gate is not satisfied; no provisional heatmap belongs in the poster.

### F4 — GFL/GFM rescue matrix

**Question:** Is the H4 mechanism specific to one technology assignment, or does a GFL/GFM substitution rescue stability?

**Panel design**

- 4 x 4 matrix.
- Rows encode the technology assignment at buses `(30, 33)` in the order `GG, GM, MG, MM`.
- Columns encode the assignment at buses `(35, 37)` in the same order.
- Each cell contains the dominant `alpha`; optional secondary text contains `g*`.
- Use the same stable/boundary/unstable semantic colors as F1.
- Add a separate marker for cells whose run is numerically invalid or failed a contract check; an invalid cell is not “stable.”

**Required data**

```text
assignment_id, row_assignment_30_33, column_assignment_35_37,
technology_by_bus, alpha, dominant_lambda, boundary_g_star,
status, MAC_to_H4, eigenvalue_gap, mode_family_id, solver_backend,
equilibrium_hash, model_hash, run_id
```

**Acceptance criteria**

- All 16 cells use the same network, loading, operating point convention, and continuation policy.
- GFL/GFM labels are explicit at the bus level; do not encode them only through color.
- A failed or missing cell is rendered as missing/invalid, never as an inferred value.

**Current readiness:** design only. Requires a dedicated Julia GFL/GFM lattice run after the F3/M1 gate decision.

## 3. Evidence strip

The strip is supporting evidence, not a fifth headline result. It contains compact, independently checkable cards:

1. **Python--Julia parity:** exact GFL11 H4 comparison with both `alpha`, dominant frequency, state count, and absolute/relative differences.
2. **TDS checkpoint:** only after a valid time-domain campaign has a declared initialization, event, sampling, and pass/fail contract. Include the triple/H4/repaired cases only if the run is complete and traceable.
3. **Reproducibility token:** short `run_id`, `model_hash`, `equilibrium_hash`, and software/backend labels. Link the detailed tables from the report, not from ad-hoc desktop files.

The current exact Python--Julia GFL11 reconciliation is eligible for the first card. The latest M1A result is diagnostic evidence and must not be presented as TDS validation.

## 4. Shared data contract

Every figure-producing run must write one self-contained run directory under:

```text
reports/poster/ias2026/research/bnd_h4_mechanism/results/<run_id>/
```

The directory must contain:

```text
config/       frozen parameters, grid, tolerances, software versions
inputs/       immutable or hashed input references
raw/          full spectra, eigenvectors, operators, and solver diagnostics
derived/      figure-ready tables only
figures/      PNG/PDF/SVG generated from this run
claims/       machine-readable claim and gate status
report/       RUN_SUMMARY.md and validation tables
```

Minimum provenance fields for every row or file:

```text
run_id, git_commit, dirty_worktree, model_hash, equilibrium_hash,
solver_backend, precision, tolerance_profile, timestamp_utc
```

For modal or reduced-order work, retain the complete eigenvalue/eigenvector set, modal ordering/tracking metadata, Schur complement residuals, complement condition number, and root-solver diagnostics. A CSV containing only `alpha` is not sufficient for a claim-bearing run.

### Optional mode-identity row for F3/F4

F3 and F4 should additionally emit the following row whenever a parameter sweep or technology substitution can change the dominant mode:

```text
MAC_to_H4, eigenvalue_gap, mode_family_id
```

These fields are interpretive diagnostics, not new gates. They distinguish “the H4 mode moved left” from “the H4 mode lost dominance and another mode took over.” `MAC_to_H4` must use the frozen H4 reference eigenvector, `eigenvalue_gap` must report the separation from the next competing tracked mode, and `mode_family_id` must remain stable under continuation unless the run explicitly records a mode-family transition.

No run may write temporary files, archives, copied source trees, or ad-hoc outputs at repository root. Temporary material belongs in a run-local `tmp/` directory or an OS temporary directory and must be removed before the run is marked complete. Poster assets belong only in the run's `figures/` directory until explicitly promoted.

## 5. Readiness matrix

| Figure | Design | Data | Scientific gate | Next permitted action |
|---|---:|---:|---|---|
| F1 Stability Cliff | ready | partial/existing | exact 16-case table | collect or reconcile the 16 spectra |
| F2 Local vs collective | ready | baseline available | baseline contract pass | regenerate figure from physical local factors |
| F3 Feedback map | ready | blocked | M1 strict pass | repair base/flagship matched load blocks |
| F4 GFL/GFM matrix | ready | not collected | common Julia lattice contract | run only after the gate decision |
| Evidence strip | ready | Python--Julia partial | TDS card requires valid TDS | add only traceable evidence |

## 6. Freeze rules for the next run

Before execution, the run must declare:

- the figure(s) it is intended to populate;
- the exact input hashes and branch commit;
- the grid, tolerances, and pass/fail thresholds;
- the output directory and cleanup policy;
- the claim status expected if a gate fails.

If a gate fails, preserve the diagnostics in the run directory, mark the affected figure `BLOCKED`, and stop downstream experiments. Do not create a provisional figure from a failed reduction, and do not package the repository into a zip as part of the run.
