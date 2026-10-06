# Q0 — Frozen ExpN/ExpP baseline reproduction

Status: **PASS_BASELINE**. Exact ExpN and ExpP candidate hashes were verified before model calls. No input candidate, frozen model, Project, or Manifest was modified.

- Rebuilt the ExpP mixed PD architecture and trimmed the declared P/Q shares. Maximum trim residual: `5.395296582857102e-12`; maximum component P/Q contract errors: `1.1368683772161603e-15` / `1.4813844595451542e-15` pu.
- Recomputed the full physical spectrum: analytic alpha `-0.0500000011050941` s⁻¹, PD alpha `-0.050000000138500054` s⁻¹, 101 physical poles.
- Recomputed the ExpP bus-38 conditional algebraic root: `1.1243444770645108` MW versus incumbent `1.124344477653483` MW.
- Rechecked the nominal frequency-zero robustness witness: beta upper `3.746493388818478e-11` < required `1.6991206999182038e-6`.
- Independently simulated the sustained 100 MW bus-16 load step for 60 s using PD. Peak frequency deviation `39.33340125578227` Hz; peak sampled RoCoF `1.072208443371192` Hz/s; settling `NaN` s. Q/P loads remain separate at buses 31 and 39.
- Runtime records are in `Q0_RESULTS.json`; event TDS used one independent PD build and one 60 s simulation.
