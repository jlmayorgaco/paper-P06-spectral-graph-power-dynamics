# Experiment G — analytic robust SG to GFL co-design

**Status: PARTIAL.** Candidate hash: 2047437a10d650b251aea5a563dba828c300e072c6913078e508e13e9d6f1b5c. Design did not import or call PowerDynamics.

## Frozen system and endpoint gates

The program discovered buses 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, reconstructed 5402.761089978776 MW of initial SG electrical injection, and preserved initialized impedance-load values including buses 31 and 39. The all-SG spectrum comparison status is **PASS** (Hausdorff residual 2.0463630789890885e-12); no-load KCL residual is 7.184918763469982e-14.

The all-GFL endpoint was classified **DEFECTIVE_JORDAN** across three deterministic gain patterns. The global-angle gauge residuals are 7.062482803481528e-19, 8.134680800347634e-19, 1.3148999464528228e-18; the physical zero remains when gains reach their frozen maxima. The scalar authority is therefore treated as a local quotient-branch seed.

## Design results

The exhaustive linearized surrogate returned **SURROGATE_GLOBAL** with support 39 and retained 32.52506590423567 MW. It is not the full-closure optimum.

The best deterministic exact-spectrum coordinate candidate used order authority_per_MW_order, retains 3821.771246228776 MW, converts 1580.98984375 MW, and has spectral abscissa -0.13793268780355464 s^-1. Its conservative modal robustified abscissa is -0.05000107709644368 s^-1, against required -0.05 s^-1.

The numerical full-block small-gain check at beta=1.6991206999182038e-6 is **true** with beta-star 0.0003700545429037017, peak resolvent 2702.304347227612, and margin 0.9954084587461466. The modal-residue separation is more conservative and is reported separately.

## Robustness, RoCoF, and transient

Retained inertia is 71395.24629937788 MVA s against requirement 6000.0 MVA s. Transient status is **EVALUATED**; unit RoCoF gain is 0.0008859066137572665 Hz/s per MW and unit frequency gain is 0.003993882658579204 Hz per MW. Direct/modal relative errors are 2.194059571812654e-9 and 3.701944542101513e-8.

The local perturbation screen status is **PASS_ON_TESTED_RETENTION_AND_RELOCATION_DIRECTIONS** over 156 deterministic directions. Full-closure KKT/SOSC is **NOT_CERTIFIED_FULL_CLOSURE_KKT_SOSC**; neither a strict local optimum nor a full-closure global optimum is claimed.

## Blinded ExpE comparison

Only after writing and hashing the pre-blind candidate did the runner read ExpE's provisional result. ExpE retained SG MW parsed as 1.1892118821653361, anchor bus parsed as NaN; the comparison and source hash are in TABLE_G14_blinded_ExpE_comparison.csv. ExpG's candidate was not changed after this read.

## Remaining independent gates

Post-freeze PowerDynamics spectral validation and nonlinear TDS are complete: candidate α=-0.13795151380393075 s^-1; both validation gates report PASS. The candidate remains unchanged and full-closure global optimality is not certified.
