# CDW68 deviation log

Every entry uses the system clock and cites the commit it follows. Preregistration commit: `6f188531` (2026-09-13T09:19). No gate, threshold, candidate, policy or model parameter has been changed.

## 2026-09-13T10:00:47-0500 (system clock read immediately before writing; after `d86037e6`): execution fixes and post-hoc labels

### Execution fixes (no definition changed)

1. **Model B converter P/Q check.**
   - ANDES 2.0 does not populate `PV.q.v` after the power flow, so the H17-style check `pq_error` read 0 and reported a spurious 2.8 pu mismatch.
   - The check was replaced by the forced-PQ consequence it guards: the converter-bus voltages of the nominal network equal the internal power flow (max 3.4e-7 pu on the full V68 case).
   - The forced-PQ targets themselves were verified (`q0` equals target at every converted bus).
2. **Model B JSON writing.** A numpy boolean in the mode records was not serializable (the first census launch wrote no file). A numpy-aware serializer was added.
3. **Model B parallelism.**
   - `multiprocessing` spawn fails inside the ANDES venv (WinError 87 in `spawn_main`, `OpenProcess` on the parent pid), and a second launch stalled.
   - The stalled processes of this campaign were stopped. Their command line contains `xtool-andes-gfl`; no foreign process was touched.
   - All ANDES work now runs as N independent single-process shards (`run68.py shard`, `launch_shards.ps1`). Each writes the same per-task JSON, so the results are unaffected.
4. **ANDES model registration.**
   - SG68D/S/M are registered at runtime, and code is generated into the private folder `home/pycode_cdw68`, seeded by a copy of `home/pycode`.
   - The shared isolated ANDES copy used by the hardening campaign is not modified.

### Post-hoc analyses (labelled exploratory; they cannot change any preregistered verdict)

5. **EM-tracked branch sensitivity (R10–R12, Model A).**
   - The preregistered critical eigenvalue of the fully converted target is the real PSS-washout pole (about −0.067 s⁻¹, gap 2e-4 to its neighbour). Its branch effects are about 1e-5 to 1e-4 s⁻¹.
   - The same perturbed matrices are therefore also projected on the rightmost EM-band mode, and its finite tracked change is recorded. Reported as "post hoc, EM-tracked".
6. **TDS holdout (R14).**
   - The preregistered selection yields no case: no level-D or level-T reversal exists, and no nested pair meets the level-D conditions needed for a negative control.
   - An exploratory TDS takes instead the two EM-tracked nested pairs with the largest same-sign marginals at the R15 reference policy.
   - It uses the same BDF-plus-network-solve scheme as `ibr_cycles.nonlinear.tds`. That tool's observer and guards read IEEE-39 device fields and cannot be applied to the 68-bus devices.
   - Reported as "post hoc".

### Interpretation notes recorded before the branch and TDS results were read

7. **Materiality of global marginals.**
   - In every REAL, NOGOV and SP33 portfolio of both models, α⊥ is a slow real pole: the PSS washout (1/T_W = 0.0667 s⁻¹), the converter Q/V leak (0.05 s⁻¹) or the common-frequency mode.
   - No global marginal reaches τ = 0.01. Levels A–D are therefore empty by construction, and the fixed-ranking FULL stratum has no material preference.

## 2026-09-13T10:31:03-0500 (system clock; after `d86037e6`, before the R8–R14 checkpoint commit): TDS timing and figure fixes

8. **Post-hoc TDS time axis.**
   - The exploratory TDS applies the reactor over 0–10 s and fits the decay over 11–29 s. The preregistered protocol is 1–11 s and 12–30 s.
   - Every run starts at the equilibrium, so the two protocols differ only by a 1 s shift of the time axis. The numbers were not recomputed.
   - The rerun at about 10:25 only added the stored traces for F10. Its decay changes agree with the six digits printed for the 10:04 run (whose file was overwritten before any commit; 10:54 amendment).
9. **R14 integrator guards.** Entry 6 already records that the generic voltage/speed guards of `tds.simulate` were not applied. Bus voltages were not stored, so guard compliance is not claimed. The tracked-mode modal amplitude stays below 5e-5 (F10).
10. **Figures.** Layout-only fixes to F6, F8, F10 and F12 (overlapping titles, axis offsets, truncated labels). No plotted value changed.

## 2026-09-13T10:38:10-0500 (system clock; after `49309c1e`): correction of entry 7, R15 execution fix

11. **Correction of entry 7.** Entry 7 said that no global marginal reaches τ in REAL, NOGOV and SP33. That is true for REAL and SP33 of both models: the largest level-A |Δα⊥| is 8.7e-3 s⁻¹ (A REAL), 4.7e-3 (SP33) and 2.3e-5 (B REAL). It is **false for NOGOV**.
    - NOGOV has three global marginals at or beyond τ: −0.0110, −0.0106 and −0.0100 s⁻¹.
    - All three replace G12, at P68_02 and P68_05, and all are stabilizing.
    - They sit on the 0.0025 Hz common-frequency mode that appears once governors are removed. They are level B (same mode, MAC ≥ 0.9999) but not level C, because the mode is outside the EM band.
    - All have the same sign, so no reversal follows, and the NOGOV verdict (0/6 at every level) is unchanged.
    - Entry 7 is kept as written. This entry supersedes it for NOGOV.
12. **R15 script.** `ranking_subset()` failed on an empty `raw/R15L` store while the ranking draws were still running. It now returns an empty table; no definition changed.

## 2026-09-13T10:54:53-0500 (system clock; after the R19 decision): paper and report production

13. **Paper edits (R23, CASE B).** Done only after all gates and the CASE decision were computed.
    - New title, abstract, contributions C1–C4, a new IEEE-68 section, and revised discussion, limitations and conclusion.
    - All IEEE-68 numbers are macros generated by `code/R23_paper_numbers.py`.
    - To keep the length at 11 pages (the previous version had 10), two IEEE-39 side results were shortened: the modal-mixing subsection now has one sentence, and the topology subsection lost its figure (old Fig. 8). Both remain complete in the unchanged supplement. No IEEE-39 number or claim changed.
    - The stacked-float page was fixed by setting one figure to 88 % width.
14. **Report conversion (R22).** The CDW converter `experiments/cdw/md2tex.py` is used unchanged, through `code/R22_report_tex.py`, which only replaces the title block and adds unicode glyph mappings. Tables that were indented inside list items were un-indented in the Markdown, because the converter does not detect them.
15. **Runtime note.** Another session's multi-worker job shared the 22-thread machine during this campaign (many foreign Python processes were visible at 10:27). This campaign used at most 16 workers and never touched a foreign process. Wall times are therefore not benchmark timings.
