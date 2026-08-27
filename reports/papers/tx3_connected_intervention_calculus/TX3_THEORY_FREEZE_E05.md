# TX3 theory freeze for E05 mechanism--consequence mapping

**Freeze ID:** `TX3-TF-E05-1.0`  
**Parent freezes:** `TX3-TF-0.2` and `TX3-TF-C1-C2-1.0`.

E05 preserves the C1/C2 action coordinates, fully re-equilibrated descriptor operator, primary logdet functional, frequency band, and direct finite Mobius convention. The frequency-resolved externality density is the signed Mobius sum of `log|det T(j omega)|` at the registered coalition vertices. Its normalized log-frequency integral must recover the E03 direct externality within the registered numerical comparison tolerance.

Dynamic modes are obtained from the exact algebraic Schur complement `A = M^{-1}(f_x-f_y g_y^{-1}g_x)`. Positive-imaginary OP00 modes in 0.1--30 Hz define immutable reference families. Families are mapped across operating points and action vertices by one-to-one Hungarian assignment maximizing scale-invariant left/right biorthogonal MAC. No `min(zeta)` mode switching is used.

For each tracked modal scalar `q`, its finite consequence externality uses the same Mobius operator as C1: `Delta_S q = sum_U (-1)^(|S|-|U|) q(1_U)`. The primary scalar is damping ratio; real part, frequency, log pole-resolvent factor, joint endpoint changes, and system spectral abscissa are secondary. A large negative logdet externality is not defined as harmful. Harm requires a negative damping-ratio interaction or positive real-part interaction on the locked family.

E05 tests `C3a`, not the original cycle-localization C3: at least one objectively selected E03 coalition has a reproducible material tracked-mode consequence under the frozen holdout gate. No quiver cycle, SCC, surgery target, correction, causality, TDS consequence, or EMT transfer is claimed.
