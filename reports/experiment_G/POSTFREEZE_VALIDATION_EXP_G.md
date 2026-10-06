# ExpG post-freeze PowerDynamics validation

Candidate SHA-256: `2047437a10d650b251aea5a563dba828c300e072c6913078e508e13e9d6f1b5c`. The validator checked this hash before loading PowerDynamics. No candidate field or controller gain was changed.

Overall post-freeze validation status: **PASS**. All-SG reference: EVALUATED, spectral abscissa -0.09806540933624897 s^-1. Candidate equilibrium: EVALUATED, residual infinity norm 1.0463605494025893e-12, 117 finite poles with 2 numerical gauge poles; candidate spectral abscissa -0.13795151380393075 s^-1, strict-margin pass true.

The analytical reduced model has spectral abscissa -0.13793268780355464 s^-1. The PowerDynamics model includes its stock inverter filter, current, DC-link and PLL states; its full finite spectrum is reported separately rather than treated as an identical state realization. The alpha difference is a cross-model comparison, not a pole-by-pole identity claim.

Nonlinear TDS status: PASS. Pulses at bus 16 use fractions 0.0005 and 0.001 of the frozen load setpoint. See TABLE_G17_powerdynamics_nonlinear_tds.csv and the trajectory table.

Relative response-scaling errors between the two load pulses are 5.2380183230615884e-5 for retained-SG COI frequency and 9.347696812600503e-6 for voltage.

The blinded ExpE comparison is in TABLE_G14. Its `source` label indicates anchor bus 38; ExpG retains 3821.771246228776 MW and fully converts bus 38. This comparison table is source-derived; the frozen candidate file remains untouched.

The PowerDynamics initializer reported `InternalLinearSolveFailed` while its residual remained within the documented acceptance tolerance; the measured candidate fixed-point residual was 1.0463605494025893e-12, and `fixed_point=true`.

No retuning or post-validation candidate edits were made.
