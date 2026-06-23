# Prior-Art Audit and Certificate Pivot

Date: 2026-06-08

Status: preliminary but action-guiding. This memo separates what is already
covered by nearby literature from the remaining plausible original direction.

## 1. Prior-Art Findings

Claude's strategic critique is mostly correct. The field is more crowded than
the current manuscript framing implies.

### C1: Damping-region / critical-mode geometry

Overlap is high. The proportional or homogeneous damping physics is already
present in power-system spectral work. A 2026 arXiv paper on frequency shaping
for oscillation damping in weakly connected power networks explicitly studies
damping ratio and decay-rate rules for inter-area oscillations:

https://doi.org/10.48550/arXiv.2601.19665

Therefore C1 should be framed as a region visualization and planning coordinate,
not as a new damping formula.

### C2: Diagonal screen plus coupling-triggered correction

This remains the most defensible direction, but it has adjacent prior art.
Pan et al. propose a modal perturbation method for non-proportionally damped
systems and use a coupling index alpha to quantify non-proportionality:

https://doi.org/10.3390/app10010341

The power-system gap is not simply "non-proportional damping exists." The gap is
whether one can provide an a priori diagnostic for the damping-margin error of a
cheap diagonal screen in IBR grids.

### C3: SG-to-IBR / GFM placement ranking

Overlap is serious. Liyanage et al. rank grid-forming inverter placement using
frequency and damping-ratio indices, including validation on a modified IEEE
39-bus system:

https://doi.org/10.1109/OJIES.2025.3538480

C3 should not be claimed as "ranking inverter placement by damping" in general.
It can remain as an application only if the ranking is produced by the proposed
cheap diagnostic estimator and is evaluated by regret against a full eigensolve.

### C4: Damping-ratio sensitivity

The QEP sensitivity formula is not new as a general tool. Setareh and
Ghazizadeh derive closed-form QEP damping-ratio sensitivities and rank remedial
actions for synchronous-generator models:

https://doi.org/10.1016/j.compeleceng.2025.110501

The differentiable piece is narrower: apply damping-ratio sensitivity to line
reinforcement/topological changes in IBR-rich grids and show disagreement with
frequency-only line ranking.

### A priori error estimation

A priori error estimation exists in grid-forming converter model reduction, but
for simulation model reduction rather than damping-margin screening:

https://doi.org/10.1049/cp.2019.0067

This supports Claude's proposed pivot: do not compete on "another estimator";
compete on "when is the cheap screen trustworthy?"

## 2. Candidate Original Core

The strongest remaining research question is:

> Can we certify, before solving the full QEP, whether the diagonal modal damping
> screen is reliable for an IBR grid with non-proportional damping?

In modal coordinates:

```text
x_ddot + Gamma x_dot + Lambda x = 0
Gamma = Delta + E
Delta = diag(delta_k)
E = offdiag(Gamma)
```

The diagonal screen ignores `E`:

```text
zeta_diag = min_k delta_k / (2 sqrt(nu_k)).
```

The exact QEP is:

```text
P(s) = s^2 I + s(Delta + E) + Lambda.
```

A useful contribution would be an a priori bound:

```text
|zeta_full - zeta_diag| <= f(E, modal gaps, diagonal poles)
```

with a trigger:

```text
if f <= tolerance:
    use diagonal screen
else:
    solve reduced QEP
```

This would be stronger and more original than simply proposing another damping
estimator.

## 3. First Numerical Audit of Candidate Bounds

Implemented in:

```text
temp/audit_diagonal_error_certificate.py
temp/diagonal_error_certificate_report.json
```

Three bound families were tested on synthetic 39-bus-like and degenerate cases.

### 3.1 Bauer-Fike companion bound

Result: mathematically safe when finite, but too loose for paper value.

Example:

- weak IBR: finite in 74% of cases, covers the error, but median tightness is
  about 1.9e6.
- mild/strong IBR: usually infinite.

Verdict: not usable as the central contribution.

### 3.2 Row/Gershgorin local bound

Result: too conservative. It usually declares "not certified" because modal
disks overlap or the discriminant fails.

Verdict: useful only as a fail-safe trigger, not as a publishable sharp bound.

### 3.3 QEP-specific second-order indicator

Because `E_kk = 0`, the first-order pole shift from off-diagonal modal damping
is zero. A second-order indicator is:

```text
|Delta s_k| <= |s_k|^2 / |2s_k + delta_k|
               sum_{l != k} |E_kl|^2 / |p_l(s_k)|
```

where:

```text
p_l(s_k) = s_k^2 + delta_l s_k + nu_l.
```

Translated to damping-ratio error:

```text
Delta zeta_k <= |Delta s_k| / (|s_k| - |Delta s_k|)
```

Numerical behavior:

- It produced no false-safe cases in the tested ensembles for tolerances
  0.0025, 0.005, 0.01, 0.02, and 0.05.
- For weak IBR, at tolerance 0.01 it certified 96% of cases with max certified
  actual absolute error about 0.00114.
- For rho = 0.20, at tolerance 0.02 it certified 45.6% of cases with no false
  safe decisions.
- For mild, strong, and near-degenerate cases, it mostly refused to certify,
  which is conservative and consistent with "solve reduced QEP."

Limitation:

- It is still loose and has weak correlation with actual error.
- It is an a priori conservative diagnostic, not yet a tight theorem.

## 4. Honest Verdict

The pivot is valid, but not yet proven enough to make the paper title a
"certificate" paper.

The current evidence supports this weaker but defensible claim:

> A QEP-specific second-order off-diagonal damping indicator can act as a
> conservative a priori trigger: it certifies low-error diagonal screening in
> weakly non-proportional IBR regimes and refuses certification in coupled or
> near-degenerate regimes.

The stronger claim still needs work:

> A rigorous, nontrivial, tight a priori error bound for the diagonal damping
> screen.

To reach that stronger claim, the next mathematical step is a small-gain or
Schur-complement proof around each diagonal pole, using the second-order term
and a remainder bound. The next empirical step is to run the same certificate on
ANDES IEEE 39-bus linearizations.

## 5. Recommended Paper Repositioning

Do not headline C1 or C3 as novelty.

Reposition as:

1. C2 core: coupling-diagnostic diagonal screen with a conservative a priori
   QEP error trigger.
2. C4 core/application: damping-aware line reinforcement for IBR grids, compared
   explicitly against frequency-only ranking.
3. C1 support: damping-region visualization.
4. C3 support: planning use case, differentiated by cheap estimator and regret
   against full eigensolve.

Suggested revised title direction:

> A Coupling-Diagnostic Spectral Screen for Damping-Margin Planning in
> Inverter-Dominated Power Grids

## 6. Follow-up Audit: Closed-Form Second-Order Correction

Claude later proposed a stronger pivot: use the second-order perturbation term
not only as a conservative trigger, but as an actual closed-form correction to
the diagonal damping-margin estimate.

This was independently checked in:

```text
temp/audit_second_order_qep_correction.py
temp/second_order_qep_correction_report.json
```

For the modal QEP

```text
P(s) = s^2 I + s(Delta + E) + Lambda,
```

with `Delta = diag(delta_k)` and `E = offdiag(Gamma)`, the diagonal pole for
mode `c` is `s0`, where:

```text
p_c(s0) = s0^2 + delta_c s0 + nu_c = 0.
```

Eliminating the other modal coordinates gives the second-order correction:

```text
Delta s_c =
    (s0^2 / (2 s0 + delta_c))
    sum_{l != c} E_cl E_lc / (s0^2 + delta_l s0 + nu_l).
```

The corrected pole is:

```text
s_c^(2) = s0 + Delta s_c.
```

This is standard perturbation logic for a QEP, but the power-system use is
specific: estimate how non-proportional inverter damping moves the damping
margin away from the diagonal inertial-Laplacian prediction.

### Independent surrogate results

The correction was tested on 300 random cases per family. Error is relative
error against the full QEP/eigensolve margin.

```text
weak_ibr:
  diagonal median / p95:      0.038% / 0.881%
  second-order median / p95:  0.007% / 0.163%
  QEP r=6 median / p95:       0.024% / 0.821%
  QEP r=20 median / p95:      0.004% / 0.093%

mild_ibr:
  diagonal median / p95:      2.53% / 33.21%
  second-order median / p95:  0.65% / 34.57%
  QEP r=6 median / p95:       1.51% / 20.30%
  QEP r=20 median / p95:      0.22% / 2.53%

strong_ibr:
  diagonal median / p95:      5.16% / 66.10%
  second-order median / p95:  1.46% / 82.86%
  QEP r=6 median / p95:       3.90% / 47.22%
  QEP r=20 median / p95:      1.24% / 24.75%

near_degenerate:
  diagonal median / p95:      307.59% / 568.37%
  second-order median / p95:  74.35% / 520.51%
  QEP r=6 median / p95:       159.98% / 328.95%
  QEP r=20 median / p95:      9.77% / 36.69%
```

Interpretation:

- The second-order correction is real and directionally correct.
- It beats the diagonal estimate in most cases:
  - weak IBR: 92.7% of cases,
  - mild IBR: 86.0%,
  - strong IBR: 76.0%,
  - near-degenerate: 91.3%.
- It often beats a small reduced QEP (`r=6`) in median error.
- It does **not** reliably beat a larger reduced QEP (`r=20`) in this
  independent run. Claude's stronger claim that second order beats `r=20`
  should be treated as single-source until code and sampling are reconciled.
- In the tail, second order can worsen p95 for strong or clustered regimes.
  Therefore it needs a validity trigger; it should not replace QEP universally.

### PLL / first-order system claim

Claude also reported a 7-state PLL verification with median improvement
`2.16% -> 0.09%` over 399 cases. After the model was provided, it was versioned
and locally reproduced in:

```text
temp/pll_7state_model.py
```

Canonical result:

```text
zeta_min (full 7-state) = 0.02905
order-0  = 0.02954  err 1.68%
order-1  = 0.02927  err 0.77%
order-2  = 0.02899  err 0.23%
```

Monte Carlo with the provided seed:

```text
399 stable cases
order-0 : median 2.16%  p95 29.2%
order-1 : median 0.93%  p95 27.1%
order-2 : median 0.09%  p95 12.9%
```

Additional audit:

- order-2 improves over order-0 in 95.5% of cases,
- order-2 improves over order-1 in 92.0% of cases,
- order-2 worsens over order-0 in 18/399 cases,
- the worst order-2 relative error is about 406%, so a trigger/fallback is still
  mandatory.

The canonical pole tracking is physically coherent. The full critical pole is
approximately `-0.034236 +/- j1.177983`; the order-0 critical pole is
`-0.034827 + j1.178430`; second order moves it to
`-0.034154 + j1.177846`, reducing the complex-pole distance from about
`7.41e-4` to `1.60e-4`.

This upgrades the PLL result from "single-source" to "locally reproduced on the
provided reduced PLL model." It still does not imply validation on ANDES
REGCA/REGCP or other vendor-grade IBR models.

### Revised verdict

The second-order correction is more promising than the earlier bound-only
pivot. A defensible new core is:

> A closed-form second-order correction to the diagonal damping-region estimate,
> used with a coupling/degeneracy trigger that falls back to a reduced QEP.

The paper should **not** claim:

> second order always beats reduced QEP.

The paper can claim, after further validation:

> second order is a cheap analytic correction that substantially improves the
> diagonal screen in weak-to-moderate non-proportional damping regimes and
> provides interpretable modal-coupling terms; reduced QEP remains necessary for
> clustered or strongly coupled cases.

For the reduced PLL model, the analogous claim is:

> second-order perturbation of the critical pole around a memory-free PLL
> approximation substantially improves the damping-margin estimate and provides
> a quantitative red-vs-control separation; the result is reproduced locally on
> the provided 7-state model but remains pending on ANDES-grade inverter models.
