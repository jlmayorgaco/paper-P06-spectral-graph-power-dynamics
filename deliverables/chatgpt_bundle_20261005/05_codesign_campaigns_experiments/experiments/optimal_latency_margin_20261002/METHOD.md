# Mathematical method and scope

## Model and fixed contract

The source IEEE-39 nonlinear DAE, generator dispatch, nodal replacement semantics, gain bounds, and PLL measurement-error delay channel are reused without change. The two input nodal designs are SHA-256-frozen in `FROZEN_PROTOCOL.json`. For the replacement-only experiment, five nodal ρ vectors interpolate along their stored path, while every Kp and Ki is held at the seed's nominal nodal vector (`FIXED_GAIN_DESIGNS_FROZEN.json`). The criterion is the existing physical spectral margin, **Re s ≤ −0.05 s⁻¹**.

For frozen y=(ρ,Kp,Ki) and physical delays τ, the full linear DDE characteristic is

\[
\Delta(s;y,\tau)=sE-A_0(y)-\sum_i A_{\tau_i}(y)e^{-s\tau_i}.
\]

The previously verified action identity factors the fixed-model descriptor design difference as \(U\Theta(y,s,\tau)V^H\), with at most 30 action columns in the 282-variable descriptor and 203 physical modes. The determinant lemma produces a small return ratio; the delayed channel has ten PLL actions. Those are **algebraic identities of this declared model**, not claims that every root of a DDE has been found in an unbounded half-plane by symbolic algebra. The supporting independent campaign is `../analytical_delay_codesign_mega_20261002/M0_ACTION_SPACE_REPRODUCTION.csv` and `Q6` validation. This experiment imports its action/contour implementation rather than repeating the campaign.

The contour oracle counts zeros in the margin-relevant right half-plane through the trace integral of the action-space determinant and an independently accumulated phase/winding count. It retains \(e^{-s\tau_i}\), with no Padé spectral surrogate. A SAFE/UNSAFE call requires both counts to agree near an integer, controlled phase step, and small quadrature diagnostic; numerical failures become INDETERMINATE. A model-derived finite contour radius is used by the source implementation. That radius and the numerical integration are not interval-certified. Critical roots are refined against the exact characteristic and followed as delay varies. The first **observed** crossing has safe/unsafe full-contour brackets ≤0.1 ms; local root refinement supplies the finer displayed estimate.

At a simple root with left/right null vectors, the exact NEP local derivative is

\[
\frac{\partial\lambda}{\partial p}
=-\frac{w^H\Delta_p(\lambda)v}{w^H\Delta_s(\lambda)v},
\qquad
\partial_{\tau_i}e^{-s\tau_i}=-s e^{-s\tau_i}.
\]

The ten selected \(\partial\operatorname{Re}\lambda/\partial\tau_i\) values were checked against centered finite differences. The Kp, Ki, and ρ sensitivities in `T03` are **finite-difference estimates of a relinearized full model**; they are not mislabeled analytic gradients. At a simple local boundary, the implicit prediction \(d\tau_c/dK=-\nabla_K\operatorname{Re}\lambda/(\partial_\tau\operatorname{Re}\lambda)\) would be a tuning direction, but the gate stopped its optimization stage.

## Numerical finding and falsification

The stored seed's limiting family is at 5.26249 Hz and the stored high-ρ zero-delay-tuned design's at 5.64012 Hz. Their physical right-vector MAC is 0.00182, so treating both as one mode is inappropriate. Holding gains fixed along the ρ path keeps MAC 0.99997 between endpoints and changes the crossing by only +0.00343 ms. Thus the attractive 1.99316 ms stored-design gap is mostly associated with gain differences, not the tested replacement-only variation. The gain swap is merely an admissible comparison and establishes no optimal τcrit*(ρ).

No assertion is made about monotonicity outside 87.5–88.45514%, arbitrary nodal ρ directions, heterogeneous delays, larger τ intervals, global gain optima, or nonlinear delayed disturbance survival. Any global infimum defining τcrit would additionally require exclusion of earlier crossings throughout the continuum, rather than sampled contours and a tracked branch. The practical outcome is the preregistered action-space fallback plus a fast-mode diagnostic, not the intended replacement-frontier poster.
