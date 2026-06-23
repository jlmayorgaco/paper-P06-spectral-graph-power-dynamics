# Phase Gate: Inertia Assignment over Controller-Aware Bridge

This gate asks whether inertia-placement diagnostics survive when evaluated over the
controller-aware Schur NEP bridge.  The script perturbs only documented inertia
parameters present in the case files.  In the current IEEE39 IBR cases, REGF1 does
not expose a clear `M`, `H`, or `Tj` virtual-inertia parameter, so GFM virtual-inertia
assignment remains blocked unless a physical parameter mapping is supplied.

- Commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`
- Random seed: `20260608`
- eps: `0.01`

## Case Results

### no_pss

- Critical full mode: zeta=0.154474, freq=1.37049 Hz, family=I_network, pi_c=0.6804.
- Critical bridge mode: zeta=0.154474, freq=1.37049 Hz, family=I_network, pi_c=0.6804.
- Gate-1 bridge-vs-ANDES: **PASS**; sign_rate=1.000, top3_overlap=1.000, finite_signal=10/10.
- Gate-1 local NEP sensitivity: sign_rate_vs_full_matched=1.000, top3_overlap_vs_full_global=1.000, median_rel_error_vs_full_matched=0.0254.
- Gate-2 substitutability: **NETWORK_MODE_SCOPE**. Critical mode is network-family; structural substitutability can be tested on retained network modes.
- Gate-3 inertia sign reversal: **SIGN_REVERSAL_OBSERVED_ON_FULL_ANDES**; negative_full=5/10.
- Virtual-inertia validation: **BLOCKED_NO_EXPLICIT_GFM_VIRTUAL_INERTIA_COLUMN**.

### base

- Critical full mode: zeta=0.154429, freq=1.37056 Hz, family=I_network, pi_c=0.6526.
- Critical bridge mode: zeta=0.154429, freq=1.37056 Hz, family=I_network, pi_c=0.6526.
- Gate-1 bridge-vs-ANDES: **PASS**; sign_rate=1.000, top3_overlap=1.000, finite_signal=10/10.
- Gate-1 local NEP sensitivity: sign_rate_vs_full_matched=1.000, top3_overlap_vs_full_global=1.000, median_rel_error_vs_full_matched=0.0254.
- Gate-2 substitutability: **NETWORK_MODE_SCOPE**. Critical mode is network-family; structural substitutability can be tested on retained network modes.
- Gate-3 inertia sign reversal: **SIGN_REVERSAL_OBSERVED_ON_FULL_ANDES**; negative_full=5/10.
- Virtual-inertia validation: **BLOCKED_NO_EXPLICIT_GFM_VIRTUAL_INERTIA_COLUMN**.

### mix60

- Critical full mode: zeta=0.123566, freq=13.7417 Hz, family=II_control, pi_c=0.8532.
- Critical bridge mode: zeta=0.123566, freq=13.7417 Hz, family=II_control, pi_c=0.8532.
- Gate-1 bridge-vs-ANDES: **WEAK_SIGNAL**; sign_rate=1.000, top3_overlap=1.000, finite_signal=0/4.
- Gate-1 local NEP sensitivity: sign_rate_vs_full_matched=1.000, top3_overlap_vs_full_global=1.000, median_rel_error_vs_full_matched=0.00017.
- Gate-2 substitutability: **NOT_APPLICABLE_TO_CONTROL_CRITICAL**. Critical bridge mode is control-family. Projecting the network-inertia shunt residual onto this control-family mode would not test the original red-inertia classifier.
- Gate-3 inertia sign reversal: **NO_SIGN_REVERSAL_OBSERVED**; negative_full=0/4.
- Virtual-inertia validation: **BLOCKED_NO_EXPLICIT_GFM_VIRTUAL_INERTIA_COLUMN**.

## Verdict

- Overall status: **CONSERVATIVE_SCOPE**
- Abstract implication: Do not claim validated virtual-inertia assignment over the controller-aware bridge. Report GENROU.M sensitivity checks and keep GFM virtual inertia as blocked pending a physical parameter mapping.

## What this gate did not test

- It did not validate a real GFM virtual-inertia knob because the tested REGF1 sheets do not expose one.
- It did not prove that the red-inertia substitutability classifier applies to control-family critical modes.
- It did not replace the registered estimator/planning validation over the NEP bridge.