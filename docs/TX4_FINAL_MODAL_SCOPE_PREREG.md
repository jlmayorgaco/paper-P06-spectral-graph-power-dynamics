# TX4 final modal-scope audit — preregistration

## Frozen provenance and scope

This preregistration is committed before the new numerical audit. It starts
from `2cbb860eb4ba053edee99c1a040c8024f3f6828e` on a new branch and preserves
all earlier TX4 outputs as read-only inputs.

## Sign convention

The collective normalized closure is

```text
C_H(s) = I + Q_H(s)
det(C_H)=0  =>  -1 is in sigma(Q_H).
```

For `H={i} union R`, the contextual return is

```text
R_i|R(s) = Q_iR(s) (I+Q_RR(s))^(-1) Q_Ri(s)
det(I+Q_H)=det(I+Q_RR) det(I-R_i|R)
det(I-R_i|R)=0  =>  +1 is in sigma(R_i|R).
```

The two unit conditions are equivalent under regularity of the proper-subset
factor, but their eigenvalue signs are opposite because the two operators are
different.

## Modal targets

The targeted reduced-predictor band is fixed at

```text
0.3 Hz <= |Im(lambda)|/(2*pi) <= 1.5 Hz.
```

For each frozen V9 portfolio, after equilibrium and transverse reduction:

```text
alpha_EM     = max Re(lambda) in that band
alpha_global = max Re(lambda) over all finite transverse eigenvalues.
```

Structural neutral modes are excluded exactly as in the authoritative
transverse classifier: the rotational/frequency center is removed by the
transverse quotient and residual numerical zero modes with `abs(lambda) <=
1e-3` are not treated as physical critical modes. No other eigenvalue is
discarded.

The global critical-mode taxonomy is frozen as:

```text
APERIODIC_REAL   : |Im(lambda_global)|/(2*pi) <= 1e-6 Hz
SLOW_OUTSIDE_BAND: 1e-6 < f_global < 0.3 Hz
TARGET_EM_BAND   : 0.3 <= f_global <= 1.5 Hz
FAST_OSCILLATORY : f_global > 1.5 Hz
```

The slow false-safe audit includes every global false-safe with
`1e-6 < f_global < 0.3 Hz`; it is selected by frequency/mode class, not by
alpha magnitude. A clean oscillatory slow case is eligible for optional TDS
only if its equilibrium, eigenvalue, and mode diagnostics are finite and
physically interpretable. At most three are selected by increasing
cardinality and then lexicographic portfolio label.

## Physical local-factor test

At each point of the frozen scalar `g` sweep, construct the pre-normalized
physical action matrix `M(s)=D(s)K(s)` in the fixed rectangular port basis.
For each H4 device record

```text
L_i(s) = I + M_ii(s), sigma_min(L_i), abs(det(L_i)).
```

Also record `sigma_min(I+Q_S)` for H4 and all 15 proper subsets at the H4
boundary frequency. The collective interpretation passes only if all four
physical local factors remain nonsingular over the recorded boundary
neighborhood while the H4 collective factor approaches zero.

## Frozen classification gates

- V4 sign audit: all four `Q_H` values must approach `-1` and all four
  contextual returns must approach `+1`, with finite Schur residuals.
- Physical local-factor gate: PASS only if no `I+M_ii` factor becomes
  singular before the collective H4 factor in the boundary neighborhood.
- V9 global task: retain the frozen 395/512 result as a separate global task.
- V9 targeted task: report confusion metrics against `alpha_EM >= 0`; this is
  never called a global safety classifier.
- P7: remove from the poster if its unstable mode is outside the validated
  approximately 0.62-Hz collective mechanism; otherwise retain with explicit
  modal scope.
- No result is upgraded to universal or probabilistic generalization.
