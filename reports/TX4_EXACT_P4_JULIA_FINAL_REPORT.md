# TX4 Exact P4 / GFL11 Julia Reproduction

## Verdict

The exact frozen TX4 P4 flagship closes as Case A for the requested custom same-model Python-vs-Julia reproduction.

The independent Julia implementation reproduces the frozen Python V4 census with the exact 11-state GFL, exact P4 SG policy, matched AC equilibrium, 16-case verdict agreement, and critical voltage-mode MAC of `0.9999999999999988` for H4.

This is an IAS poster-credibility upgrade for the stated same-model spectral claim. It is not an EMT, PowerDynamics SimpleGFLDC, or TPWRS-equivalence result.

## Frozen setup

- H4: `{30,33,35,37}`.
- P4: `(g,k,t,h)=(0.03625,1.425,1.5,1.0)`.
- GFL11 state order: `theta_pll, x_pll, p_filt, q_filt, x_p, x_q, i_d, i_q, x_id, x_iq, x_v`.
- Replaced devices use matched P/Q initialization, so all 16 cases share the frozen AC point.
- The prior GFL10 harness remains separate and is not relabeled as exact.

## Primary gates

| Gate | Result |
|---|---|
| G1 exact GFL11 | PASS; device-level derivatives and injections match exactly at the fixed parity point, max difference `0.0`. |
| G2 H4 state count | PASS; `86=4*11+6*7`. |
| G3 16-case verdict census | PASS; 16/16 verdicts identical. |
| G4 H4 alpha | PASS; Python `0.1270064680506376`, Julia `0.12700646832229723`, difference `2.72e-10 s^-1`. |
| G5 H4 frequency | PASS; Python `0.6222796695029233 Hz`, Julia `0.6222796695187043 Hz`, difference `1.58e-11 Hz`. |
| G6 H4 voltage-mode MAC | PASS; `0.9999999999999988`. |
| G7 equilibrium | PASS; maximum voltage-magnitude difference `2.29e-14`, maximum angle difference `2.30e-08 rad`. |
| G8 hidden mismatch audit | PASS; equation and P4 policy audits are explicit and device/all-SG parity files are present. |

All 15 proper subsets are stable in both codes. H4 is unstable in both codes with the same positive transverse critical mode.

## State-count correction

The prior independent Julia H4 result had 82 states because it used GFL10. The exact reproduction has 86 states because four replaced converters each add `x_v`.

## Evidence

- `results/TX4_EXACT_P4_V4_CROSSCODE.csv`: complete 16-case census.
- `results/TX4_EXACT_P4_ALLSG_PARITY.csv`: all-SG parity.
- `results/TX4_GFL11_DEVICE_PARITY.csv`: device-level equation parity.
- `results/TX4_EXACT_P4_MODE_MATCH.csv`: H4 modal match.
- `results/TX4_EXACT_P4_STATE_COUNTS.csv`: GFL10/GFL11 state-count distinction.
- `figures/TX4_F1_ALPHA_CROSSCODE.png`, `TX4_F2_FREQUENCY_CROSSCODE.png`, `TX4_F3_H4_VOLTAGE_MODE.png`, `TX4_F4_STATE_COUNTS.png`.

## Limitations

The result is a numerical reproduction of the declared reduced custom IEEE-39 DAE and its index-one spectral linearization. It does not validate EMT behavior, a second independent dynamic model, controller saturation, protection, faults, or TPWRS submission readiness. The exact allowed claim is therefore limited to cross-language reproduction of the frozen P4/GFL11 same-model spectral result.

## Reproducibility

Branch: `research/tx4-exact-p4-julia-reproduction`.

Parent: `069fa21204aa7c6f45fb44eb17e304ab99685a2b`.

No push was performed.
