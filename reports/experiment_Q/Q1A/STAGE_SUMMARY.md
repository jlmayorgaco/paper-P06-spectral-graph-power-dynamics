# Q1A — Frequency-output definitions and metric audit

Status: **METRIC_DEPENDENT**. Frozen ExpP candidate only; no design variable was changed.

- Primary metric frozen for Q: maximum absolute deviation among every present SG rotor-speed output and every present local GFL PLL-frequency output on generator buses 30–39. Rotor speed uses `60(ω−1)` Hz; the installed PLL state is rad/s, converted by `Δω/(2π)`.
- Independent bus diagnostic: derivative of unwrapped bus-voltage phase, using a cubic 11-sample Savitzky–Golay local polynomial at 100 Hz. Window support 0.1 s; measured -3 dB derivative-response bandwidth `11.70691475` Hz. It injects no plant signal.
- PD traces cover sustained 1, 10, and 100 MW steps and a 100 MW, 0.1 s bus-16 pulse. Per-family peaks, RoCoF, limit decisions, and bus/time locations are tabulated.
- Frequency metric classification is `METRIC_DEPENDENT` because frequency families produce at least one different inherited-limit decision. The primary rotor/PLL metric is fixed by the specification; bus frequency remains an independent diagnostic.
- Wall time was not instrumented. The stage contains four independent PowerDynamics TDS cases; each trace uses the stored 0.01 s sample interval.
