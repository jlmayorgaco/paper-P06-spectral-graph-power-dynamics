# E32 — nonlinear phasor-domain DAE time-domain validation: **PASS**

This is **not** an EMT study. It is a **nonlinear phasor-domain DAE simulation**:
the network stays algebraic, the devices keep their differential equations, and
nothing is linearized. It asks whether the modes the eigenanalysis reports are
the modes the nonlinear system actually exhibits.

Seven configurations × four disturbances, **25 runs completed**, 3 not applicable
(D4 steps a converter at bus 30, which the base case, the single-replacement case
and the keep-30 restoration do not have).

Disturbances, identical magnitude across every configuration:

| | disturbance |
|---|---|
| D1 | +2 % active load at bus 20 |
| D2 | +2 % active load at bus 29 |
| D3 | **+0.2 %** mechanical power at the machine on bus 31 |
| D4 | +2 % active-power reference at the converter on bus 30 |

D3 carries its own magnitude, declared and identical across cases: a 2 %
mechanical step on a single machine drives the speed spread out of the
small-signal window before a ringdown can be fitted, so it would not have been a
small-signal test at all.

The step is applied as an initial-condition mismatch rather than a discontinuity
inside the right-hand side — the formulation E25 arrived at after the stiff
solver stalled on the switched version.

## Result

| quantity | value |
|---|---|
| matrix-pencil frequency error, median | **0.00066 Hz** |
| matrix-pencil frequency error, worst | **0.0029 Hz** |
| matrix-pencil damping error, median | 0.0053 |
| matrix-pencil damping error, worst | 0.030 |
| **damping sign agreement, matrix pencil** | **25 of 25** |
| unstable runs whose ringdown grows | **4 of 4** |
| stable runs whose ringdown decays | **21 of 21** |

The flagship is predicted unstable at `+0.1447` and 0.5747 Hz. The nonlinear
ringdown returns `+0.1347` to `+0.1452` at 0.5744 to 0.5757 Hz across all four
disturbances. The repairs are predicted at `−0.0498` (RC) and `−0.0500` (RB) and
measured at `−0.0498` to `−0.0805`. The condenser is predicted at `−0.2894` and
measured at `−0.2812` to `−0.2897`.

**RC reconstruction check.** RC was rebuilt from the frozen E21 record, whose
detail line omits `tau_p`; it was taken unchanged and the reconstruction verified
before use: spectral abscissa `−0.0498` against the recorded `−0.0500`.

## Two estimator caveats, reported rather than smoothed over

**The crude log-envelope regression disagrees on one run of 25** — the flagship
under D3, where it returns `−0.0007` while the matrix pencil returns `+0.1452`
against a prediction of `+0.1447`. A single exponential cannot represent a signal
carrying several modes at a peak amplitude of `2.3e-5`; the growing inter-area
component has not yet dominated the decaying ones inside the 12 s window. The
primary damping statistic is the matrix pencil, the estimator E25 validated for
exactly this reason, and it agrees on all 25 runs. The envelope number is
reported beside it rather than dropped.

**The FFT peak is a single-peak statistic** and picks a different band mode on
two runs (RC and RB under D3, 1.367 Hz against a 0.559 Hz prediction), which is
why its worst error is 0.81 Hz against a median of 0.040 Hz. The pencil resolves
the modes; the FFT is kept as the independent cross-check it was asked to be, and
its failure mode is stated.

## What this establishes

The linear eigenanalysis on which the entire Track-A argument rests is confirmed
by nonlinear integration of the same DAE, for the flagship instability, for both
converter repairs, for the condenser and for machine restoration, under four
independent disturbances. Growth and decay go the way the linear model says, at
the frequency the linear model says, to under 0.003 Hz.

This is an internal consistency check between two code paths in the same project.
It is not an independent implementation; that is E31, which is **G2 MIXED**.

## Files

`E32_TDS_summary.csv` (+ `.parquet`), `E32_TDS_frequency_validation.csv`,
`E32_flagship_vs_repair_TDS.png`, `E32_figure_source.csv`, `manifest.json`.
