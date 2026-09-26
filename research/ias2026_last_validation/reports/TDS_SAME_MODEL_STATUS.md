# Same-model Julia TDS

status: `NUMERICALLY_VERIFIED_WITH_OBSERVABILITY_LIMIT`

The fresh Julia fixed-step nonlinear DAE traces cover the base, proper
30+33+35 subset, repaired three-target case, four-target blocker, and a
critical-mode-seeded blocker. The common 2% load pulse and algebraic Newton
solve use the exact frozen custom SG/AVR/PSS/GFL equations. Maximum retained
algebraic residual is below `1e-9`.

The base/proper/repaired traces are consistent with their stable transverse
spectra. The common short pulse did not visibly excite the RHP mode in the
four-target trace, and the mode-seeded trace did not produce a robust
independent nonlinear growth estimate over the retained horizon. Therefore the
blocker is not promoted as a nonlinear-TDS proof; its primary evidence is the
cross-code reconciled spectrum and the Julia collective mechanism audit.

Raw traces: `raw/tds_same_model/`.
