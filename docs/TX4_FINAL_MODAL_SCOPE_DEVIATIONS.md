# TX4 final modal-scope audit — deviations

Initial status: `NONE_RECORDED`.

Any deviation from the frozen scope must record the affected rule, reason,
exact action, stop/continue decision, and claim impact. Numerical results and
negative outcomes will not be edited to satisfy a desired case.

## Deviation D1 - legacy V9 table was not same-policy data

- Affected rule: the preregistered separation between the frozen global V9
  task and the targeted 0.3--1.5 Hz EM task, with the requirement that any
  retained global comparison be on the same frozen P4 policy.
- Reason: audit of the discovery scripts showed that the archived
  `TX4_V9_BLIND_VS_FULL.csv` was generated through the FC10/default matched-q
  path, while the blind predictor and this audit use the frozen TX4 P4 policy
  (`g=0.03625`, `k=1.425`, `t=1.5`). The policy difference changes portfolio
  classifications, so the archived 395/512 and 116 false-safe values are not
  a same-policy benchmark.
- Exact action: recompute all 512 V9 full transverse spectra under the P4
  policy after preregistration. Use the resulting 332 stable and 180 unstable
  global cases, 110 global false-safe cases, and the separately reported
  511/512 EM-band result as the authoritative numbers for this audit. Preserve
  the legacy CSV as a frozen historical input and label it policy-confounded.
- Stop/continue: continue the finite modal-scope audit; do not rerun any
  broader campaign or modify the frozen predictor.
- Claim impact: no claim may cite 395/512 or 116 as a same-policy result. The
  corrected report must disclose the conflict and must not imply that the
  targeted EM classifier is a global safety classifier.

## Deviation D1 - legacy V9 table was not same-policy data

- Affected rule: the preregistered separation between the frozen global V9
  task and the targeted 0.3--1.5 Hz EM task, with the requirement that any
  retained global comparison be on the same frozen P4 policy.
- Reason: audit of the discovery scripts showed that the archived
  `TX4_V9_BLIND_VS_FULL.csv` was generated through the FC10/default matched-q
  path, while the blind predictor and this audit use the frozen TX4 P4 policy
  (`g=0.03625`, `k=1.425`, `t=1.5`). The policy difference changes seven or
  more portfolio classifications, so the archived 395/512 and 116 false-safe
  values are not a same-policy benchmark.
- Exact action: recompute all 512 V9 full transverse spectra under the P4
  policy after preregistration. Use the resulting 180 unstable and 332 stable
  global cases, 110 global false-safe cases, and the separately reported
  511/512 EM-band result as the authoritative numbers for this audit. Preserve
  the legacy CSV as a frozen historical input and label it policy-confounded.
- Stop/continue: continue the finite modal-scope audit; do not rerun any
  broader campaign or modify the frozen predictor.
- Claim impact: no claim may cite 395/512 or 116 as a same-policy result. The
  corrected report must disclose the conflict and must not imply that the
  targeted EM classifier is a global stability classifier.
