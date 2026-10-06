# Q2 feasibility reference — all SG

This is an upper-bound reference, not an optimized Q design. It uses the frozen ExpN model, controller bounds and the exact +100 MW bus-16 sustained event with Q unchanged.

- Complete physical spectrum: alpha=-0.09806540933624713 s⁻¹ across 113 physical poles.
- Robustness: beta requirement 1.6991206999182038e-6; omega=0 pointwise upper witness 0.007691096114689081. The Bounded Real CARE/LMI result is BOUNDED_REAL_CERTIFIED, beta lower bound 1.6991206999182038e-6, CARE relative residual 1.8635950338134124e-6, minimum P eigenvalue 0.0024072583421757094, maximum CARE/LMI residual eigenvalue -0.0840422912292011, and closed-loop abscissa -0.04806540772094596. It uses Float64 arithmetic without directed rounding, so it is numerically certified, not formally validated.
- Analytic step model: F_inf=0.03755683117282334 Hz, F_peak=0.07580838269071084 Hz, R_peak=0.12532351438704908 Hz/s. Independent PD TDS over 60 s: local SG frequency peak 0.07536489665023405 Hz, SG RoCoF 0.12383507394416611 Hz/s, bus-frequency diagnostic peak 0.07464181551095384 Hz.
- Trim residual 2.6645352591003757e-13, P/Q maximum error 0.0 pu. The result JSON records 45.664 s for TDS and 171.968 s total for the measured run; an earlier stage note contained a stale TDS timing.
- Reference feasible under declared constraints: true. The all-SG fallback candidate was frozen before PD/TDS validation; no optimization or globality claim is made.

The comparison trace at `TABLE_Q2_all_sg_linear_vs_PD.csv` has 6001 samples. The maximum absolute difference between the max local rotor-frequency envelopes is `0.0009813 Hz`; the analytical and PD peaks are `0.0758081 Hz` and `0.0753649 Hz`. This validates this event/output for the all-SG reference only.
