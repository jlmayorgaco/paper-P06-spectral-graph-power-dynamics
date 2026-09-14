# TX4 contextual-return deviations

This file is initialized before new campaign numerics.

No deviations are authorized at preregistration time. Any deviation must be
appended with date/time, frozen rule affected, reason, exact action, whether
the gate stopped, and claim impact. A computational convenience may not change
the model, disturbance, threshold, coordinate, or inclusion rule.

Initial status: `NONE_RECORDED`.

## Pre-valid numerical implementation erratum

Before accepting any campaign result, an exploratory invocation of the new
script used `I+R` in the numerical Schur residual while the frozen TX4 block
partition requires `I-R`. The exploratory output was invalidated and is not
used as a result. The preregistration and implementation were corrected to
the exact identity in the TX4 plan; the valid campaign is rerun from the
corrected code. This changes no model, data, grid, threshold, or claim.

The initial implementation list duplicated triple 30+33+35 as both a
predecessor and a control. Before accepting the valid run, the list was made
nine-case and unique by replacing that duplicate control with singleton 33.
No numerical outcome or threshold was used to choose it.
