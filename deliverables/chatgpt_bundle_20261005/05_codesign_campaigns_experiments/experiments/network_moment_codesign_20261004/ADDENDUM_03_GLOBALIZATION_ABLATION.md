# Globalized corrector and fair ablations

The complete fourteen-pattern multimode run is preserved in TABLE06--09.
Uniform+1ms still failed (alpha=+0.05507/s) after60 iterations: a5% step bound
alone did not prevent alternating active modes and rejected linear predictions.
This is not evidence of infeasibility. No failed result is overwritten.

New outputs live in `globalized/`. Use the same bounds, exact moment equality,
eleven tracked modes,60-iteration budget and5% step limit, adding a standard
merit backtracking check. Merit is half squared relative gain effort plus1000
times positive violation of the -0.0600001/s catalog guard. Backtrack over
1,1/2,...,1/1024 until Armijo1e-4 holds (roundoff allowance1e-8). Track every
trial's roots. Termination update norm1e-6 and guard error1e-7; these are numerical
termination tolerances, not a relaxation of the independent -0.05/s pass gate.
Archive nonconvergence/failure, continue to other patterns.

Repeat collective co-design for all fourteen original patterns. For the fixed
uniform+1ms headline also run the identical optimizer without the moment
equality and the closed-form moment-only minimum. This measures the constraint's
cost and the need for modal protection. No superiority or maximum replacement
claim follows from this adaptive algorithm-development experiment.
The original full-spectrum/five-event headline is applied to the resulting
uniform+1ms collective candidate only if it passes its catalog gate.
