# Preregistered graph/GSP joint-design experiment

Frozen before the new model computations. Date: 2026-10-03. No commits/pushes;
all new work stays here and previous frozen results remain read-only.

Question: can graph representations reduce controller design dimension without
losing feasible SG-to-GFL replacement, under identical local search conditions?
This is a bounded local predictor/corrector comparison, not global optimization.

Baseline: rho_i=0.875, Kp_i=0.9*(2*pi*5), Ki_i=(2*pi*5)^2/4, uniform pure PLL
error delay 40 ms. These are fixed before evaluating this new design. If baseline
is not feasible, report the failed gate; do not change it to force a positive result.

Variables: all ten rho_i are free in EVERY family. Only the representation of
log Kp/log Ki differs. Frozen gain limits remain 0.25--4 times nominal.
Families: common; low physical-graph modes q=1,2,3,4,5; full nodal gains;
high physical-graph q=3; seeded random q=3 (seed 20261003); dynamic-port q=3.
Every reduced basis contains the constant vector. Physical graph is the frozen
lossless synchronizing Kron graph on candidate buses 30--39. Dynamic-port basis
uses the three lowest eigenvectors of Re(D_H) at the preregistered 5 Hz and 40 ms,
orthogonalized against the constant vector. No choice uses trial-design labels.

Primary comparison: one identical linearized joint-design subproblem from the
same baseline, with exact local DDE root gradients and discrete nonlinear delayed
event tangents. Trust bounds are |delta rho_i|<=0.015 and |delta log K_i|<=0.08.
Maximize initialized GFL dispatch. Among objective ties within 1e-7 fractional
dispatch, minimize squared nodal gain movement. Correct each proposal through
fixed fractions 1, 1/2, 1/4, 1/8; evaluate all five events for each tested design.
Record all failures. A fully nodal space contains every reduced space: a reduced
family cannot beat its exact global optimum; a possible benefit is efficiency
or better nonlinear behavior for this finite local search.

Security: required spectral abscissa <=-0.05/s; all five frozen events
(8,-100), (16,+100), (16,-100), (29,+100), (29,-100) MW; post-event horizon 60 s;
frequency/RoCoF windows 0.5 s, both limits 0.5 in Hz and Hz/s; voltage [0.9,1.1];
normalized SG actuator slack >=0.002. No current-limiter safety claim.
Simulating the event at t=0 with an equilibrium history is a time translation
of the frozen t=1 event and its 61-s record. Pre-event voltage history is retained.

Tangents: second-order SDIRK method of steps, h=0.02 s (tau/h integer), with exact
derivatives of its discrete delayed trajectory. Accepted points use independent
adaptive Rodas5P method of steps, tolerance 1e-9 and maximum step 0.01 s. Tangent
versus adaptive baseline errors must be <=0.002 Hz, <=0.002 Hz/s and <=2e-5
actuator fraction; otherwise refine to h=0.01. Validate one seeded parameter
direction by centered differences of the discrete method before optimization.
If any guard/branch is nonsmooth, its derivative is not promoted to a smooth proof.

Full numerical DDE contour count is required before calling a returned design
spectrally feasible. Exact exponential equations are used; floating-point counts
are not interval certificates. Report event feasibility separately if spectral
coverage is unresolved. Source-model/parametric-Jacobian parity gate: <=1e-8 relative.

Graph expansion test, independent of optimization: at baseline, frequencies
0.1,0.5,1,2,4,5,8,12 Hz and delays 0,20,40 ms. Test diagonal graph-modal impedance
and Neumann corrections of orders 0,1,2,3,5,8 against the full impedance inverse.
Convergence bound is permitted only when ||D^-1 E||_2<1; otherwise report failure
or descriptive error. These are walks in the modal interaction graph, not
Markov probabilities or demonstrated physical cycle causality.

A GSP benefit requires all-event feasible replacement within 1 MW of the full
nodal result with fewer variables. Superiority over random/high-frequency spaces
is reported separately; no cherry-picking among graph modes. Methodological
novelty, global optimality, general heterogeneous-delay effects and a complete
maximum-replacement frontier are outside what this finite comparison can prove.
