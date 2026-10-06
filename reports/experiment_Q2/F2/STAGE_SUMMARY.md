# F2 — ExpN grid-frequency modal residues

The full physical quotient spectrum of the frozen ExpN point was expanded at buses 30–39. Finite pole residues use q=(C r)(lᴴB)/(lᴴr); the direct algebraic feedthrough D is separately retained. The phase jump from the step produces a distributional component and is not represented by the finite-pole residue sum.

- Smooth filtered-free grid frequency peak (state response plus D): 21.254654352603954 Hz at bus 34, t=60.0 s after the step.
- Smooth finite-state RoCoF peak: 12334.913423686541 Hz/s at bus 39, t=0.0 s; ideal unfiltered RoCoF also has an impulse because phase jumps.
- Steady offset magnitude: 23.10444870318589 Hz at bus 34.
- Active finite mode: -0.0500000011050941 + 0.0im. Residue condition number: 17758.87465876337; six-time modal reconstruction error: 1.5287326959878555e-11 Hz; zero-time RoCoF reconstruction error: 7.075868779793382e-10 Hz/s.
- Contribution rows contain per-mode residue and contributions evaluated at each metric's independently determined maximum.

The ExpN rightmost mode dominates the frequency peak and steady offset, but it is not the dominant RoCoF mode. At the RoCoF maximum, the pair at `-5622.35 ± 875.02j s⁻¹` contributes about `-6367.14 Hz/s` per pole; other fast modes partially cancel it. At the frequency peak, the rightmost pole contributes `+38.51 Hz` and the next slow real pole contributes `-17.41 Hz`. For the steady offset those contributions are `+40.53 Hz` and `-17.58 Hz`.

For this same +100 MW event the measured algebraic bus-phase jump is `0.0012662112 rad` at bus 38, so ideal unfiltered electrical-frequency RoCoF is distributional/unbounded. The sampled 100 Hz finite-difference RoCoF is `2.01815 Hz/s`; the Savitzky–Golay sensitivity values for 0.04, 0.10, 0.20, and 0.40 s windows are `1.39992`, `0.80353`, `0.78338`, and `0.74261 Hz/s`. All finite windows tested exceed the frozen 0.5 Hz/s project limit at this ExpN point. The frequency peak is `21.25465 Hz` and steady offset `23.10445 Hz`. These are analytical full-model values; no 100 MW PD bus-voltage measurement was run after the new estimator was introduced.
