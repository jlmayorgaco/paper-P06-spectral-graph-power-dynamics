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
   - The rerun at about 10:25 only added the stored traces for F10. Its decay changes are bit-identical to the 10:04 run.
9. **R14 integrator guards.** Entry 6 already records that the generic voltage/speed guards of `tds.simulate` were not applied. Bus voltages were not stored, so guard compliance is not claimed. The tracked-mode modal amplitude stays below 5e-5 (F10).
10. **Figures.** Layout-only fixes to F6, F8, F10 and F12 (overlapping titles, axis offsets, truncated labels). No plotted value changed.
