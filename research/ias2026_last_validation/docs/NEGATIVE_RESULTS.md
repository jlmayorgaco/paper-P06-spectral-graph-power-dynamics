# IAS2026 negative-results register

This file prevents negative or incomplete evidence from being silently turned
into a positive Vancouver claim.

1. The official PowerDynamics P2 and SimpleGFLDC P5 models remain negative
   relative to the frozen custom blocker: their fresh V4 spectra are stable
   after gauge-aware filtering.
2. P5 rating equivalence is not established. The corrected scope is matched
   scheduled P/Q on the common 100-MVA system base; retired-machine MVA values
   are reported only as provenance metadata.
3. The same-model nonlinear TDS traces solve with small algebraic residuals,
   but the common short load pulse does not visibly excite the frozen custom
   RHP mode. The RHP verdict is therefore carried by the reconciled spectrum,
   not by an overstated nonlinear-trace claim.
4. The second-model PLL/current bandwidth search found no mixed stable/unstable
   policy region on its frozen representative grid. A mixed >=12-point
   holdout was not fabricated after this negative result.
5. No common physical uncertainty envelope, robust radius, or new scaling sweep
   was executed.
6. No paper or poster was rewritten from this validation bundle.

These are results, not missing documentation. The raw files and reports keep
the evidence auditable.
