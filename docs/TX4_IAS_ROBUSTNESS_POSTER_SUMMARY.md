# IAS2026 Vancouver Poster Summary - TX4 Robustness

## Headline

In the frozen IEEE-39 matched-policy reduced DAE, exact H4 evaluations cover 20,590 declared conditions. The bounded engineering screen reports H4 presence of 57.8% by QMC and 58.1% by MC; exact-H4 minimality is 7.3% and 7.1%, respectively, with proper-subset rows outside the exact calibration skeleton explicitly labeled surrogate-tier.

## Safe claim

The result supports a benchmark- and policy-conditioned robustness screen, not a physical robust radius or universal IEEE-39 claim. Raw alpha_all values are preserved; blocker flags use transverse alpha_EM to remove numerical neutral-mode contamination.

## Evidence

- 329,440 master portfolio rows; all 16 portfolios per condition.
- 20,590 exact H4 rows.
- 6,816 exact all-portfolio calibration rows.
- 302,460 calibrated surrogate proper-subset rows.
- Julia cross-code: exact frozen nominal 16-case parity retained; random 16-case spot-check not executed.

## Do not claim

No physical uncertainty envelope, robust radius, EMT/current-limit/DC-link/protection/hardware result, fresh nonlinear same-model TDS, or random Julia parity.
