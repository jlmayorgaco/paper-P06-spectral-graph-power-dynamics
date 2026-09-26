# TX4 Robustness Statistics Audit — Deviations Log

1. The source campaign used legacy full-spectrum endpoint logic and then a
   transverse `alpha_EM` normalization. These columns are retained as audit
   inputs; corrected H0-based endpoints are newly generated.
2. The source campaign has no serialized ExtraTrees model predictions for the
   exact all-16 calibration conditions. Surrogate validation will therefore
   use a clearly labeled reconstructed, grouped held-out evaluation rather
   than treating in-sample refitting as independent evidence. This cannot
   qualify a surrogate for primary minimality if the strict gate is not fully
   satisfied.
3. No U005/H005 engineering threshold was found in the existing preregistration
   or TX4 code at audit start. No threshold is inferred from the name “005”.
4. Exact fallback computation, if triggered, uses the existing condition
   coordinates from the source master table and does not regenerate samples.

