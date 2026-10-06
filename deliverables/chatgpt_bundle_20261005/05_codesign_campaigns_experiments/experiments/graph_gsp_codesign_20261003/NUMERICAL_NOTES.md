# Numerical log

The first delayed-tangent invocation failed before producing results because a
nested observer reused the name of the zero-error history vector as a state
index. Renamed both bindings. This was an implementation correction; the frozen
experiment definition, baseline, graph families and acceptance limits did not change.

The preregistered h=0.02 tangent mesh exceeded the 2e-5 actuator-fraction tolerance in four events (up to 4.82e-5). Preserved those outputs in grad/baseline_h002 and refined to h=0.01 as specified, before optimization.

The h=0.01 mesh retained a 2.44e-5 actuator error at bus 29 / -100 MW. The near-linear mesh error exposed a discontinuity-handling bug: the stiffly accurate endpoint stage at t=tau used the post-jump delayed detector when integrating the preceding interval. Corrected this single endpoint to its left limit; the next stage still uses post-event history. Old meshes are preserved. Repeat both mesh and derivative gates before any optimization.

Before nonlinear candidate evaluation, SLSQP failed to resolve the declared minimum-gain-norm tie break for three families. Replaced only that convex QP numerical solver with CLARABEL at fixed constraints/objective; preserved all initial outputs in proposals_initial_solver. This avoids accepting arbitrary LP extreme points in the families whose tie-break solver failed.

The first parallel candidate batch stopped before evaluating any design because independently included Julia modules had distinct DesignContext types. Passed the matching ReducedDAE module explicitly to the spectral wrapper. Physical matrices and algorithm are unchanged; no failed numerical candidate was discarded.

The full and low_q1 first proposals lie very near the -0.05/s contour. Adaptive trace integration hit depth 15 with local errors 6.31e-7 and 7.01e-7 (threshold5e-7); these were unresolved computations, not accepted spectra. Increased only the quadrature refinement allowance to22, retaining the same contour and acceptance thresholds. Added explicit INDETERMINATE recording and checkpoint restart.
