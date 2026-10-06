# Free-pattern gain geometry — theory before numerical screening

Question: can the exact full-network PLL return and bounded real PI gains
exclude a prescribed complex frequency for every modal pattern and every
gain vector in a declared box? Can such exclusion close a continuous contour
and support a gain-uniform design decision?

This is a semidefinite relaxation of exact physical gain equations, not a new
general SDP theorem or a replacement-maximum claim. Fixed delays remain 40 ms.
The preceding interaction campaign is immutable input.

First derive row-wise Hermitian inequalities, including a gain-parallelogram
enclosing disk. Every true eigenpair must satisfy them. A nonnegative weighted
sum that is negative definite excludes all patterns. Numerical SDP output is
only a candidate certificate until interval arithmetic verifies the full
matrix and any continuous frequency panels. Gauge and eliminated hidden
states require separate treatment. A feasible SDP is inconclusive.

Frozen numerical screen:
- Common replacement rho = .875, .90, .95, .975, .99.
- Real part = -.05; positive frequencies = .25,.5,1,2,3,4,5,6,8,10 Hz.
- Full historical gain bounds: Kp in [.25,4]*(2*pi*5),
  Ki in [.25,4]*(2*pi*5)^2/4, independently at each of ten sites.
- Secondary radius screen at the previous corrected and violating joint
  designs: relative boxes 0,.001,.005,.01,.025,.05 around all twenty gains.
  Rectangle boundary: Re in [-.05,.02], Im in [30.7,31.2], eight equal
  segments per side. This contour is a numerical diagnostic, not a proof.
- True-root checks: all six preceding endpoint designs and their tracked
  certified roots. The exclusion oracle must never exclude these roots.
- Diagonal-only network return is a declared ablation, never ground truth.

Primary falsifiers: a true root is excluded; all gain boxes admit the relaxed
SDP on relevant contours; point certificates cannot be extended continuously;
no engineering decision changes; a simpler established bound is equally useful.
Any adaptive follow-up gets a separate addendum, without changing these grids.

The independent causal audit tests the existing blind-interval idea and retains
its negative result. It does not define the primary spectral screen.

Claim statuses: EXACT_IDENTITY, PROVED_REDUCED_MODEL, NUMERICALLY_VALIDATED,
SUPPORTED_LOCAL, SUPPORTED_EMPIRICAL, NEGATIVE_RESULT, BLOCKED.
