# F3 — PLL gains and sustained power support

A centered finite difference over 30 withheld mixed designs perturbed Kp at 15 points and Ki at 15 points, one local gain per point. DC output is exact `-C A⁻¹B + D` from the fixed-support reduced model.

- Largest sampled `|dH0/dKp|`: 7.656890899819576e-14 Hz/MW per gain unit.
- Largest sampled `|dH0/dKi|`: 1.2459340800613709e-14 Hz/MW per gain unit.
- Largest sampled `|dα/dK|`: 7.593602048107098e-8 s⁻¹ per gain unit.
- Dynamic centered differences at the actually active GFL bus 38 in the corrected ExpG interior point are in `TABLE_Q2_F3_dynamic_gain_derivatives.csv`: Kp sensitivities are α 6.02×10⁻¹³, pointwise β radius −3.95×10⁻⁹, h=5 Fpeak −3.71×10⁻⁸, h=5 Rpeak −4.78×10⁻⁵, and F∞ −1.36×10⁻¹⁵; Ki sensitivities are 3.34×10⁻¹³, −4.12×10⁻¹⁰, 1.92×10⁻⁸, 1.14×10⁻⁶, and 7.73×10⁻¹⁸, respectively. The β calculation is a pointwise singular-value derivative at the frozen ExpG peak frequency, not a reoptimized or formally certified β-star.
- Classification: **PLL_DC_SUPPORT_INDEPENDENT**. β is a sampled resolvent estimate, not a formal frequency-domain certificate.
