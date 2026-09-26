# E40 — numerical stress test of the closure: **the grid, not floating point, is the error**

43 cases: 20 near-boundary and 20 generic points from the E33 deterministic
continuation, which is exactly reproducible, plus the flagship, the RC repair and
the 25 % condenser at the nominal point. No explicit inverse of `T₀` is formed
anywhere.

**Extended precision is unavailable and was declared so before the run.** numpy's
`longdouble` has a **53-bit mantissa on this platform**, identical to `float64`,
and nothing may be installed. The perturbation ensembles and the three-algorithm
comparison are the substitute; the limitation is stated, not hidden.

## Floating point is not the problem

| check | worst over 43 cases |
|---|---|
| backward residual of the `T₀` solve | **3.6e-17** |
| LU vs QR vs SVD spread in `m₄` | **4.4e-14** |
| perturbation 1e-14, worst shift in `m₄` | 2.4e-12 |
| perturbation 1e-13, worst shift in `m₄` | 3.7e-11 |
| perturbation 1e-12, worst shift in `m₄` | **2.8e-10** |
| perturbation 1e-14, worst shift in the closure eigenvalue `μ` | 2.4e-12 |
| perturbation 1e-13, worst shift in `μ` | 3.7e-11 |
| perturbation 1e-12, worst shift in `μ` | **2.8e-10** |
| `cond(T₀(jω))` | 3.0e+03 |
| smallest `m₄` audited | **0.00488** |

The response is **linear in the perturbation amplitude** across three decades,
which is what a well-conditioned problem does, and the eigenvalue moves by
essentially the same amount as its distance to −1. The smallest margin audited
exceeds the worst numerical shift by a factor of **1.7e+07**.

## The frequency grid dominates, by eight orders of magnitude

| | |
|---|---|
| worst variation across the 61-, 121- and 241-point grids | **2.3e-02** |
| worst variation in `f_port` across those grids | 1.0e-02 Hz |
| worst reduction from golden-section local minimisation | 2.6e-02 |
| worst floating-point effect from the table above | 2.8e-10 |

**The grid is about 10⁸ times more important than floating point.** That is the
honest characterisation of the frozen metric: it is a *grid-resolution-limited
upper bound* on the true band minimum, and its looseness is concentrated near the
boundary, where the margin is smallest and the minimum sharpest. Refined, the
smallest audited margin falls from 0.00488 to **0.00116**.

This bias runs **against** every claim made from `m₄`: the frozen metric
overstates the margin exactly where a small margin is being claimed. It was not
replaced, and the verdicts of E33 and E35 rest on it.

## Bug found and fixed during the run

The perturbation ensemble first computed the closure eigenvalue from the
*unperturbed* operator, so the reported eigenvalue variation was identically
zero — a metric that measured nothing. It now reads the eigenvalue off the
perturbed solve. The margin statistics were unaffected; only the eigenvalue
column was wrong, and it is corrected here.

## Files

`E40_numerical_stress.csv`, `manifest.json`.
