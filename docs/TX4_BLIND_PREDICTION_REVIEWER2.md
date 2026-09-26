# Reviewer 2 — model adequacy and nonlinear validation audit

## Question

Is the PowerDynamics-derived TX4 model realistic enough, are assumptions
credible, is H4 physically meaningful, and does TDS support the interpretation?

## Assessment

The critical mode is a 0.622-Hz oscillatory inter-area mode. The modal audit
shows dominant SG rotor/EMF participation with nonzero GFL control
participation; the change is therefore best described as a network-mediated
electromechanical interaction, not a pure local controller pole. Equilibrium
residuals, independent finite-difference Jacobian checks, NumPy/SciPy
eigensolver agreement, and contextual local-factor regularity show no evidence
of a numerical/index pathology in the reported witness.

The scope remains limited: the custom ten-state GFL model is phasor-domain,
the campaign has no EMT portfolio validation, controller/current limits are
not a general hardware certification, and E31 ANDES does not reproduce the
custom GFL boundary. The frozen +2% bus-20, 0.2-s nonlinear phasor pulse gives
proper stable, H4 unstable, and retuned H4 stable, with 3/3 growth-sign
agreement. That supports the local interpretation but is not a broad transient
stability claim.

Scores (0--10): novelty 7; correctness 8; model adequacy 5; evidence 7;
usefulness 7; IAS impact 8; TPWRS readiness 4.

## Recommendation

Use the result for IAS communication with explicit model labels. For TPWRS,
add an independently reproduced custom-GFL model, limits, more networks, and
pre-registered holdouts. Preserve the negative V9 and ANDES results.
