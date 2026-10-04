# Real-gain sectors after the polygon screen

The polygon/disk relaxation excluded 30/50 sampled points for the full gain
box; the remaining points were inconclusive. Small gain boxes around the
previous pair designs also become inconclusive near the left contour edge.
These are numerical diagnostics, not stability or infeasibility proofs.

Strengthen the same necessary-condition relaxation by retaining normalized
real gain variables through their actions on each detector. Add interval
sector products, phase equalities, and four RLT cross products per site.
Lift [q;u][q;u]* and drop rank one. This is a real structured-uncertainty/IQC
relaxation, not a new general control theorem. The exact nonlinear rank-one
formulation remains equivalent to gain-box root existence on the regular chart.

Test the same 50 primary points and secondary contours for positive radii.
For every numerical witness recompute its largest Hermitian eigenvalue after
rounding nonnegative weights. Root counterexamples remain mandatory.
If useful, bracket the relative gain radius at selected frequencies by
monotone bisection; report only point-wise certificates until continuous
frequency verification and the eliminated hidden states are addressed.

New solver code is hashed before execution; no original screen is overwritten.
