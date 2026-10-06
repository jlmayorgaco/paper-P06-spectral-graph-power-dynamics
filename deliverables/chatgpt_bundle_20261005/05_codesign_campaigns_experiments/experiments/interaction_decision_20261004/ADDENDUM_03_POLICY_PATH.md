# Declared-policy replacement path, before its evaluation

The selected30/37pair is now a certified modal decision reversal. Determine
the first margin loss along its specific independent tuning policy, not an
optimum over rho allocations or all local gains.

For h in[0,.012], replace delta rho=h at both buses30 and37; prescribe each
Kp=Kp_anchor*(1+2h); solve each individual Ki so its selected pole keeps the
anchor real part. Keep all other parameters and tau fixed. Evaluate the13
coarse points h=0,.001,...,.012. Save all failures. Bisect only the first
resolved bracket crossing Re lambda=-.05 to width1e-8 in h. Verify the
two final endpoints with the same interval root method. A sign change brackets
a crossing on the numerically tracked branch; it does not certify monotonicity
throughout the entire interval, nor a global maximum.

Compute the interaction curvature at h=0 with symmetric h=1e-4 and5e-5 before
using the crossing location. Because both singleton real parts are held fixed,
the joint first derivative vanishes; I(h)=c h²+O(h³). If c>0, sqrt(b/c) is a
local predicted crossing, NOT a bound without a remainder enclosure. Compare
prediction with the executed path and report the error and loss of validity.
The corrected h=.01 design remains frozen; do not modify it using these data.
