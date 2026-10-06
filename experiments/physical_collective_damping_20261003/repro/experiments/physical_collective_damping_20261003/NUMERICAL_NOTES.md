# Numerical implementation log

1. The first spectral invocation stopped before publishing tables: the old
   40-ms root seeds alone converged to only three distinct roots for N. Added
   a deterministic frequency-seed scan from 2 to 12 Hz in 0.1-Hz increments,
   identical for all designs, and delay continuation steps <=1 ms. No physics,
   design, delay case, metric, or favorable-result criterion was changed.
2. Initial output revealed root-family jumps with 1-ms continuation and direct
   refinement between cases; those pilot tables are rejected and overwritten
   within this new, unfrozen experiment. Corrected to analytic uniform-delay
   sensitivity prediction, <=0.1-ms steps, overlap >=0.95 and adaptive step
   halving. Added the previously frozen critical roots explicitly as seeds. No
   result from the first root-following implementation is accepted. Nonlinear
   work using its generated inputs was interrupted and is restarted after repair.
