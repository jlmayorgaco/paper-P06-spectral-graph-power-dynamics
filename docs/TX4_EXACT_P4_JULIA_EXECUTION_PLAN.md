# TX4 Exact P4 / GFL11 Julia Reproduction - Execution Plan

1. Freeze the current-state record and this preregistration in a commit before running the exact case.
2. Audit the Python and Julia equations, including the `xv` state, initialization, current injection, state ordering, and P4 SG scales.
3. Run the exact Python oracle for base, H4, and all 16 V4 subsets with one fixed output root.
4. Run the exact Julia custom IEEE-39 harness for the same 16 cases, using the same frozen data and central-difference reduced Jacobian.
5. Compare verdicts, state counts, alpha, frequency, equilibrium voltages/angles, and voltage-mode shapes. Use declared Hungarian matching with frequency proximity and voltage-mode MAC.
6. Write the device-level parity and all-SG parity tables before interpreting the V4 census.
7. Run only the preregistered optional three-point `g` spot check if and only if the primary gates pass. Do not refit or select a new value.
8. Generate the four required validation figures, final report, handoff, headline JSON, and the two upload zips.
9. Render and inspect the final PDF, verify archive contents and hashes, commit the completed branch, and perform no push.

## Explicit exclusions

No SimpleGFLDC run, no EMT run, no new TDS campaign, no new theory search, no parameter tuning, and no additional model family.
