# Experiment Q2B — final summary

**Status: FAIL_OPTIMIZATION.** Best-found, PowerDynamics-validated candidate; KKT and robust interval certificates remain open.

- Retained SG: **438.446 MW** across buses 30–39; GFL dispatch **4964.315 MW** (91.885%).
- Frozen candidate SHA-256: `f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce`.
- Analytic: α=-0.050036705 s⁻¹; β observed=1.699121139343e-06 against 1.699120699918e-06; Fpeak(T=.5)=0.498425 Hz; F∞=0.491719 Hz; analytic R(T=.5)=0.126984 Hz/s.
- Robust interval: `INCOMPLETE_INTERVAL_BOUND` at 50,001 nodes. Candidate has only an observed beta margin, not a completed robustness certificate.
- KKT: primal 0.00e+00; stationarity 0.830; complementarity about 1.48e-12; LICQ rank 1; SOSC not tested.
- Independent PD: α mismatch 9.02e-12 s⁻¹; max pole mismatch 4.62e-10 s⁻¹; trim residual 4.47e-11; nonlinear bus16/100 MW Fpeak=0.4738 Hz, Rpeak=0.1259 Hz/s; both project event limits pass.
- All four tested causal RoCoF windows pass at this point, and the four out-of-sample PD events pass. This does not make it a locally optimal design.
- Support/global completeness: open. No global claim. No push or commit.

See [REPORT_EXP_Q2B.md](REPORT_EXP_Q2B.md) and the bundled CSV/figure evidence.
