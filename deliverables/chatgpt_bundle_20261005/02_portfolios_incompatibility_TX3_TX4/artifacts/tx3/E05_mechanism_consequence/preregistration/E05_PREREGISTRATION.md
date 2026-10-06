# TX3 E05 mechanism--consequence preregistration

**Freeze IDs:** TX3-E05-CANDIDATES-1.0 / TX3-E05-NUMERICS-1.0 / TX3-E05-GATE-1.0  
**Parent evidence:** C1=SUPPORTED and C2=SUPPORTED.  
**Scope:** E05 only. No ParaEMT, TDS surgery, negative control, repair, E06, or later stage is permitted.

E05 does not search for a dangerous cycle. It maps three objectively selected, reproducible E03 coalitions to externality-density spectra and BMAC-tracked small-signal modal consequences. Candidate selection is frozen from accepted E03/E04 holdout statistics before any E05 consequence is evaluated: `A7-A8` is curvature-dominated, `A3-A6` is feedback-dominated, and `A2-A7-A8` is the largest reproducible third-order coalition.

The OP00 positive-imaginary oscillatory modes in 0.1--30 Hz define reference families. After discovery only, each candidate locks the reference family having the largest median absolute damping-ratio Mobius externality among families with valid BMAC tracking on at least 12/16 discovery points. Holdout remains untouched until `E05_MODE_FAMILY_FREEZE.json` is created and committed.

A candidate is dynamically consequential only if its locked family is valid at at least 6/8 holdout points, at least 4 holdout points satisfy `|Delta_S zeta| >= 0.0025`, and the same-sign fraction is at least 0.75. The sign of the E03 logdet externality never labels an action harmful. A negative damping externality or positive real-part externality is required for a harmful interpretation.

The original cycle/SCC C3 claim remains unassessed. E05 instead tests the narrower preregistered claim `C3a`: at least one accepted finite spectral externality maps reproducibly to a material tracked-mode consequence. Negative results are retained and stop progression.
