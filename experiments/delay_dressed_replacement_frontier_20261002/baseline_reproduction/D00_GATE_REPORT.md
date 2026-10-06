# D0 baseline-parity gate

**Status: PASS.**

All three frozen design points were reconstructed from their saved rho/Kp/Ki values. Equilibrium residuals and full finite spectra were checked independently between the fixed-support ReducedDAE and PowerDynamics models. The output `TABLE_D00_FULL_SPECTRUM_PARITY.csv` contains Hungarian eigenvalue assignment errors after removing the PD rotational gauge pole.

Declared parity tolerances (absolute):
- equilibrium residual inf: 1e-08
- alpha abs error: 1e-06
- matched spectrum max abs error: 1e-05
- frequency peak abs error Hz: 0.0001
- RoCoF peak abs error Hz/s: 0.0001
- voltage extrema abs error pu: 0.0001

Observed maxima across the three design points:
- matched full-spectrum error: 5.89787e-10
- rightmost-pole error: 3.21515e-11
- frequency-peak error: 2.39896e-07 Hz
- RoCoF-peak error: 6.36953e-06 Hz/s
- Vmax error: 3.78401e-07 pu
- maximum equilibrium residual: 4.02808e-09
- five external physical holdouts complete: True

The full finite ODE spectrum, rightmost pole, and matched post-event metrics pass. The PowerDynamics validator stores one second of pre-event equilibrium before its event at t=1 s; those samples are excluded from the event-response comparison because ReducedDAE's event origin is t=0. The unfiltered full-horizon Vmax in `independent_pd/events.csv` includes this pre-event value (1.0635 pu), while the post-event maximum matches the reduced trajectory. The separate Vset=1.0 mapping remains a diagnostic variant, not the frozen event contract.

D0 compares the design event at all three stored designs and all five frozen external holdouts for the 90.047% joint candidate over matching post-event windows. See `TABLE_D00_HOLDOUT_PARITY.csv` for the five holdout errors. This is model-output parity; the joint candidate exceeds the 0.5-Hz frequency security limit in two holdouts, and actuator slack was not recomputed here. No delay-frontier experiment was run before this gate closed.
