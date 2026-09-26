# Custom-GFL reproduction in ANDES: results (Phase 8)

Date: 2026-09-11.

- Specification: `docs/20260911_GFL_REPRODUCTION_SPEC.md`, commit 3e847a4f,
  committed before any comparison.
- Preregistration: 5d0b1986, §10.
- Code: `experiments/post_cumulant_validation/gfl_repro/`.
- Outputs: `results/PCV/PCV06/`, `results/20260911_GFL_REPRODUCTION.csv`.

## Verdict

**GATE 5: PASS.** All §10 primary tolerances are met on the 32 cases (the 16
subsets of H4 at P4 and at G_S). None of the stopping rules fired.

| check | tolerance | result |
|---|---|---|
| Gate 0, device transcription (360 vectors, 18 instances) | ≤ 1e-9 (relative) | max 1.5e-14 — PASS |
| ANDES initialization residual | ≤ 1e-6 | 2.2e-13 in all 32 cases — PASS |
| bus voltages, ANDES vs internal | ≤ 1e-6 | max 5.1e-7 — PASS (see §4) |
| abs(Δα_perp) | ≤ 1e-3 s⁻¹ | max 1.26e-6 (P4 H4); every other case ≤ 1.9e-7 — PASS |
| abs(Δ critical frequency) | ≤ 1e-3 Hz | max 1.7e-7 Hz — PASS |
| RHP counts | equal 32/32 | 32/32 — PASS |
| verdicts | equal 32/32 | 32/32 — PASS |
| H, κ | equal at P4 and G_S | P4: {30+33+35+37}, κ = 4 in both. G_S: EMPTY in both — PASS |
| structural pair (two smallest abs(λ) < 1e-3, all others ≥ 1e-2) | all cases | 32/32 |
| transverse dimension (ANDES = internal n_x − 2) | all cases | 32/32 |
| secondary: nearest ANDES eigenvalue to each internal band eigenvalue | ≤ 1e-4 | max 2.9e-5 — PASS |

At P4 the complete portfolio is the only unstable subset in both tools:

| | internal | ANDES |
|---|---|---|
| α | +0.1270065 | +0.1270077 |
| frequency | 0.62228 Hz | 0.62228 Hz |
| RHP | 2 | 2 |

All 15 proper subsets are stable in both tools, with α from −0.214 to −0.144.
At G_S every subset is stable in both, and H4 has α = −0.017385 in both.

## 1. What this reproduction establishes, and what it does not

**Established.** The same device equations, solved by an independent tool,
give the same equilibrium, the same linearization and the same eigen-verdicts.
The independent parts are:
- ANDES' own network model and xlsx parser;
- the power flow;
- the initialization chain;
- symbolic (sympy) Jacobians instead of internal central differences;
- the DAE reduction;
- the eigensolver;
- a spectrum-level treatment of the rotational symmetry (the two
  structural eigenvalues removed) instead of the internal transverse
  quotient.

The P4 planning trap (H = {H4}, κ = 4) and its disappearance at G_S are
therefore **not an artefact of the internal numerical pipeline**.

**Not established.**

- **Physical adequacy of the equations.** Both tools run the same custom
  two-axis machine and the same custom GFL, by design (§4.3 of the spec).
  This is a reproduction of computation, not a validation against a
  different converter or machine model.
- **Library models.** No library GFL (for example REGCA1/REECA1) is involved,
  and no conclusion is drawn about how such models would behave.
- **Coverage.** Only P4 and G_S were run; P_inf was optional and was not run.
  No time-domain simulation was run in ANDES.

## 2. Environment and isolation evidence

- **Venv.** `.venv/xtool-andes-gfl`, from the same base CPython 3.13.14 as
  `tx3-andes`. The dependencies are pinned to the tx3-andes versions.
  `importlib.util.find_spec("ibr_cycles")` is `None`, asserted at runtime by
  `run_andes.py`.
- **ANDES.** The `git archive` of `vendor/andes` tag v2.0.0 (eda5163c),
  installed editable from `.venv/xtool-andes-gfl/src/andes`. The only change
  to that copy is:
  - the new file `andes/models/pcv_models.py`;
  - one `file_classes` line.
- **Generated code.** It lives in `.venv/xtool-andes-gfl/home/pycode`. Every
  ANDES process ran with `USERPROFILE` set to that private home and an
  explicit `pycode_path`; `run_andes.py` asserts `~ == private home`.
- **Untouched.**
  - `~/.andes`: its sha256 listing before and after the campaign is
    identical (0 differences in 303 files; `results/PCV/PCV06/env_dot_andes_*.txt`).
  - `vendor/andes`: HEAD eda5163c, 0 dirty files before and after
    (`env_vendor_*.txt`).
  - `.venv/tx3-andes`: not written; it still imports ANDES from `vendor/andes`.

## 3. Gate 0 (device level)

- **Instances.** 10 `SG2AX` instances at the P4 machine parameters
  (k = 1.425, t = 1.5), and 4 `GFL11` instances each at g = 0.03625 and
  g = 0.25.
- **Vectors.** 20 random state and terminal-voltage vectors per instance
  (seed 20260918). The states are perturbed around the internal equilibria;
  `v ~ U(0.9, 1.1)`, `a ~ U(−π, π)`.
- **Compared.** All state derivatives, and the injected `w·P` and `w·Q`
  (system base).
- **Result.**
  - maximum relative error 5.0e-15 on derivatives;
  - 1.5e-14 on `P` and 4.5e-15 on `Q`;
  - every instance passes (`results/PCV/PCV06/PCV06_gate0.csv`).

The ANDES-initialized setpoints (from ANDES' own power flow) differ from the
internal ones by at most 4.7e-7. That maximum is `pm` of the slack-bus
machine, consistent with the power-flow difference in §4. Gate 0 feeds both
sides the same setpoints, as the spec requires, so the transcription test
itself is isolated from the power flow.

## 4. Observations and deviations

- **Power-flow difference (inside tolerance, not diagnosed).** The bus
  voltages differ by at most 5.1e-7:
  - magnitudes by up to 1.7e-7;
  - angles by −1.4e-7 to −3.8e-7 rad.

  This is not a pure rotation. Both flows converge tightly (ANDES tolerance
  1e-12, initialization residual 2.2e-13), so the residual difference is at
  the data or convention level, below 1e-6. The spec requires localization only
  on a mismatch, and nothing was retuned. It is recorded as a known, bounded
  discrepancy.
- **Declared in the spec before the comparison.** A custom machine model
  (`SG2AX`), because ANDES has no fourth-order two-axis machine and
  degraded-GENROU is not an exact transcription.
- **Procedural incident (no effect on results).** The first ANDES launch
  deadlocked during ANDES' parallel code generation on Windows. The spawned
  workers re-imported the runner's top-level code, including an xlsx write.
  Recovery:
  - only that process tree (my own processes in the xtool venv) was
    terminated;
  - the partial private pycode folder was deleted;
  - code generation was re-run single-process (`andes.prepare(nomp=True)`,
    now part of `setup_xtool_env.py`);
  - the runner's side effects were moved under `if __name__ == "__main__"`.

  No comparison had been produced before the incident.
- **Not run.** P_inf (optional), and ANDES time-domain simulation.

## 5. Files

| file | content |
|---|---|
| `gfl_repro/setup_xtool_env.py` | idempotent environment builder (venv, copy, registration, private codegen) |
| `gfl_repro/pcv_models.py` | the custom ANDES models `SG2AX` and `GFL11` |
| `gfl_repro/export_internal.py` | the handoff inputs and the internal reference (tx3-analysis) |
| `gfl_repro/run_andes.py` | ANDES Gate 0 evaluation and the 32 cases (xtool-andes-gfl only) |
| `gfl_repro/compare.py` | the comparison and tolerances |
| `results/PCV/PCV06/PCV06_handoff_inputs.json` | the only file the ANDES side reads |
| `results/PCV/PCV06/PCV06_internal_reference.json` | internal Gate 0 outputs, the 32 internal cases and spectra |
| `results/PCV/PCV06/PCV06_andes_gate0.json`, `PCV06_andes_cases.json` | ANDES outputs |
| `results/PCV/PCV06/PCV06_gate0.csv`, `PCV06_gate0_summary.json` | Gate 0 comparison |
| `results/PCV/PCV06/PCV06_cases.csv`, `PCV06_comparison.json` | case comparison |
| `results/20260911_GFL_REPRODUCTION.csv` | deliverable table (one row per case) |
