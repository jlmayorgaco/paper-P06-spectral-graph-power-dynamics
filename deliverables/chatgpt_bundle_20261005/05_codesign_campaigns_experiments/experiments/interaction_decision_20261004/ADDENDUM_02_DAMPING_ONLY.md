# Damping-only versus complex-pole compensation — adaptive discovery

Frozen before evaluating this family. The previous exact complex-pole
compensation gave 0 margin violations in 230 combinations; maximum adverse
shift was only 2.91e-5 /s. Do not call that a material collective failure.

Hypothesis: preserving only the real part of an individual mode allows frequency
detuning to couple sites; preserving the entire complex pole may avoid that
specific mechanism. This is a comparison of two specified local design rules,
not a theorem that local feedback is insufficient.

Use the same anchor and target as addendum 01. At each site increase rho by
.01, set Kp/Kp_anchor=1+t for t in {-.10,-.05,-.02,.02,.05,.10}, and solve Ki
and Im(lambda) together while fixing Re(lambda)=Re(target). Keep physical gain
bounds unchanged. Solve through ten intermediate steps. Each individual action
must preserve that pole's real part to 1e-8 /s; all other spectral claims are
separate. Abort a local branch if its frequency changes by more than 2 /s.

For each t combine all 45 site pairs with the SAME t at both sites, then all
45 pairs with opposite t (positive t at lower-numbered bus). This gives 405
pair designs. Track the joint pole through ten steps. Save all successes and
failures. A false acceptance requires real-part increase taking the pole above
-.0499 /s. The nominal singleton-additive real-part prediction equals the
anchor because both individual real-part shifts are zero.

Selection: earliest |t| with a valid pair reversal; largest margin violation at
that |t| breaks ties, then lexicographic bus pair. For that pair, evaluate both
branches/mode collisions and attempt an interval root/interaction certificate.
Compare the same replacement pair with exact complex-pole compensation from
addendum 01. Any full-spectrum or nonlinear claim requires its own validation.

This is explicitly adaptive hypothesis development. No claims of independent
generalization, optimum, communication necessity, or journal novelty follow.
