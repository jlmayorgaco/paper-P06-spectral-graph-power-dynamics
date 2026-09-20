# TX4 Exact P4 / GFL11 Julia Handoff

## Final result

Case A: the exact frozen P4/GFL11 same-model spectral reproduction passes.

H4=`30+33+35+37` has `nx=86`, Python alpha `0.1270064680506376 s^-1`, Julia alpha `0.12700646832229723 s^-1`, Python frequency `0.6222796695029233 Hz`, Julia frequency `0.6222796695187043 Hz`, and voltage-mode MAC `0.9999999999999988`. All 15 proper subsets are stable in both codes, and all 16 verdicts agree.

## Scope

The closure covers only the frozen custom IEEE-39 reduced DAE, exact 11-state GFL, exact P4 SG AVR scaling, matched P/Q equilibrium, and index-one spectral reduction. It does not cover EMT, SimpleGFLDC, a second dynamic model, faults, saturation, protection, or TPWRS readiness.

## Important correction

The prior independent Julia campaign's H4 state count was 82 because it used GFL10. The exact GFL11 reproduction is 86. The old result remains separately archived and was not overwritten.

## Files

The upload bundle contains the final report, this handoff, headline JSON, claim matrix, cross-code census, modal matching, state counts, equation and policy audits, and four figures. The reproduction bundle contains source scripts, preregistration, raw tables, manifests, and deviations.

No push was performed.
