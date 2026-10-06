# Q1B–Q1G — Linear response, residues, PLL steady support and frequency bound

The inherited independent PD DAE matrix/pole identity regression passed in 69/69 ExpN cases; maximum aligned reduced-A relative error was 3.1973762568747116e-15 and maximum assigned physical pole error was 5.5648424577071397e-8 s⁻¹. These are reused from immutable ExpN evidence under the same model SHA.

The analytical +1 MW load input uses the installed IEEE-39 ZIP law at bus 16 (active and reactive fractions are entirely constant impedance; the event changes active Pset while holding Qset). The unit-step modal extrema are in TABLE_Q1C, and exact modal-residue reconstruction errors plus individual contributions of every pole at the frequency-dominant output are in TABLE_Q1D. Against independent PD trajectories over the same time window, the 1 MW and 10 MW peak errors are small; the 100 MW event departs strongly from the linear prediction.

The largest centered steady-frequency derivative over the 20 withheld gain perturbations was 1.4935078601889666e-12 Hz/(MW·gain-unit). The four dynamic modal derivative checks are -7.355798389974365e-10, 2.6291929919666317e-10, 5.057132596805252e-8, 2.2866867468514678e-8. The ExpP one-SG control perturbations test alpha, step extrema and a pointwise observed beta witness; this beta result is explicitly not a norm certificate.

The 100-point TGOV1/ZIP test classifies steady frequency as NONAFFINE; inferred K_load mean is 123.71657634497724 MW/Hz with relative range 43.92932911463566. The frequency-only continuous-knapsack output is SCREENING_ONLY. The retention allocation is not a global lower bound unless the affine structure and feasibility remainder both pass the strict test.

No optimization was run in this stage. Time histories are inherited from Q1A; no PowerDynamics call is inside an optimization routine.

Wall time was not instrumented. Work count: 69 inherited PD matrix/pole identity points, 20 steady-gain derivative points, 100 affine-stiffness points, and four dynamic derivative checks.
