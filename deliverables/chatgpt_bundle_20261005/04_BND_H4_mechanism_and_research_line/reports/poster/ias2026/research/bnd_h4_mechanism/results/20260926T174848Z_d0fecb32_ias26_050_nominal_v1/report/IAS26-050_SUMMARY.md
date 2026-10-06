# IAS26-050 — fixed intervention nominal DAE check

Status: **PASS**

Treatment was frozen before this one-point run in `reports/poster/ias2026/research/bnd_h4_mechanism/configs/IAS26-050_INTERVENTION_V1.json`: H4 all-GFL, `g=0.25`; control `g=0.03625`. No sweep, MC, reduced-root calculation, eta, or figure was run.

- Complete-spectrum `alpha_perp`: `0.127006467828` → `-0.0173846760616` s^-1 (`Delta=-0.14439114389`).
- Target-band `alpha_Omega`: `0.127006467828` → `-0.0173846760616` s^-1 (`Delta=-0.14439114389`).
- Frozen target-band Delta-alpha prediction: `-0.14439114389` s^-1; direct-DAE comparison error: `0` s^-1.
- Candidate equilibrium residuals: f `5.64599e-13`, g `7.50067e-13`; state/transverse dimensions `86/84`.
- Candidate `alpha_perp` verdict at frozen tau_dec=1e-08: **STABLE**.

Full spectra and eigenvectors for the candidate and the frozen original H4 point are retained. The direct-DAE result is not a validation of reduced poles. IAS26-060 remains prohibited until its operational distribution, correlations, limits, redispatch, seed, and 1,000 scenario IDs are frozen in its own pre-run config.
