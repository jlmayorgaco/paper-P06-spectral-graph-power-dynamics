# PD39 short follow-up: penetration versus composition and numerical audit

Date: 2026-09-14  
Branch: `research/pd39-robust-transition-confirmatory`  
Frozen discovery commit: `61554336`  
Confirmatory baseline commit: `902403cf`

## Executive result

The follow-up does not erase the 255+1 result.  All 28 six-of-eight and all
eight seven-of-eight portfolios in the frozen discovery set pass the
`m₉ >= 0.05 s⁻¹` requirement, while V8 is the only portfolio that fails it.
V8 is also truly unstable in four of the nine controller scenarios.

The new matched-composition comparison qualifies the interpretation.  The
cliff is not explained by converted MW alone: portfolios with the same
cardinality and essentially identical converted MW, IBR MVA, and remaining
inertia can differ in worst-case alpha by 0.0283 s⁻¹ for 6/8 and 0.0522 s⁻¹
for 7/8.  At the same time, the short experiment does not establish that V8's
entire cliff is a pure composition effect.  V8 remains the unique 8/8 case and
the extreme penetration/resource endpoint.  The defensible conclusion is
therefore **penetration/cardinality creates the cliff, while composition
modulates the margin within a cardinality**.

The numerical audit finds a real numerical-conditioning hazard, but not a
model-independent numerical explanation for V8.  Tight 1e-10 and 1e-12
equilibria agree; reduced-matrix and descriptor spectra agree to approximately
1e-10 or better.  A deliberately loose 1e-8 solver tolerance can produce
spurious near-zero alphas because the model is ill-conditioned, so the tight
solver setting is mandatory for interpretation.

## Frozen comparison design

The protocol is in
`docs/PD39_COMPOSITION_NUMERICAL_PREREG.md`.  It uses all eight 7/8
predecessors, V8, and six 6/8 matched pairs.  Matching used only static
metadata, not alpha or modal outcomes.  The nine existing discovery
controller scenarios are reused as frozen records; no new holdout or search
radius was introduced.

The raw outputs are:

* `results/PD39_COMPOSITION_SELECTION.csv`;
* `results/PD39_PENETRATION_COMPOSITION.csv`;
* `results/PD39_PENETRATION_COMPOSITION_PAIRS.csv`;
* `results/PD39_PENETRATION_COMPOSITION_PAIR_SUMMARY.csv`;
* `results/PD39_NUMERICAL_TOLERANCE_SWEEP.csv`;
* `results/PD39_NUMERICAL_EIGENSOLVER_AUDIT.csv`;
* `results/PD39_V8_FINITE_PERTURBATIONS.csv`.

## Penetration/cardinality evidence

The frozen discovery cardinality summary is:

| Cardinality | Portfolios | Worst alpha across 9 scenarios (s⁻¹) | Minimum m₉ (s⁻¹) | Robust cases | True-stable cases |
|---:|---:|---:|---:|---:|---:|
| 0/8 | 1 | -0.0980654093 | 0.0980654093 | 9/9 | 9/9 |
| 1/8 | 8 | -0.0979285029 | 0.0979285029 | 72/72 | 72/72 |
| 2/8 | 28 | -0.0977814965 | 0.0977814965 | 252/252 | 252/252 |
| 3/8 | 56 | -0.0976282047 | 0.0976282047 | 504/504 | 504/504 |
| 4/8 | 70 | -0.0974765740 | 0.0974765740 | 630/630 | 630/630 |
| 5/8 | 56 | -0.0979634350 | 0.0979634350 | 504/504 | 504/504 |
| 6/8 | 28 | -0.1032091539 | 0.1032091539 | 252/252 | 252/252 |
| 7/8 | 8 | -0.0858347015 | 0.0858347015 | 72/72 | 72/72 |
| 8/8 | 1 | +0.1140950119 | -0.1140950119 | 0/9 | 5/9 |

The eight 7/8 worst-case margins are:

| Missing SG | Converted MW | m₉ (s⁻¹) |
|---:|---:|---:|
| 30 | 4370 | 0.1262289835 |
| 32 | 3970 | 0.1380223155 |
| 33 | 3988 | 0.1382955669 |
| 34 | 4112 | 0.0858347015 |
| 35 | 3970 | 0.1363561762 |
| 36 | 4080 | 0.1380576591 |
| 37 | 4060 | 0.1383254195 |
| 38 | 3790 | 0.1380197164 |

V8 converts 4620 MW, corresponding to 91.89255% IBR MW penetration in the
stock system, with 407.611 MW of remaining SG dispatch in the static table.
The all-7/8 set spans 3790–4370 converted MW and still remains robust.  Within
the 6/8 set, converted-MW versus m₉ correlation is -0.2123; within the 7/8
set it is -0.3489.  These are descriptive finite-set effects, not population
claims, and they are weak compared with the abrupt 8/8 separation.

## Matched composition evidence

All selected pairs satisfy the frozen matching limits.  The strongest
within-cardinality effects are:

| Pair | Cardinality | ΔMW | ΔIBR MVA | Δremaining inertia | Maximum absolute alpha difference (s⁻¹) |
|---|---:|---:|---:|---:|---:|
| k6_p2 | 6 | 0 | 0 | 100 | 0.0283206639 |
| k6_p5 | 6 | 0 | 0 | 100 | 0.0070639042 |
| k7_p3 | 7 | 322 | 400 | 850.2 | 0.0521850149 |
| k7_p1 | 7 | 382 | 200 | 1340.0 | 0.0120665884 |
| k7_p2 | 7 | 0 | 0 | 100 | 0.0016661397 |
| k7_p4 | 7 | 20 | 0 | 210.0 | 0.0002677670 |

Four of the six 6/8 pairs have maximum |Δalpha| below 1.2×10⁻⁴.  The two
largest shifts are therefore not a universal composition offset; they are
site-specific changes that can matter despite matched aggregate ratings.
None of the selected pairs changes true-stability or 0.05-robustness class.

The figure `figures/PD39_confirmatory/F14_penetration_vs_composition.png`
(also PDF/SVG) visualizes the 6/8, 7/8, and V8 separation and selected pair
links.

## Numerical audit

The tight solver results are stable.  For V8, the nominal alpha is
`-0.00144477797` at 1e-10 and `-0.00144477803` at 1e-12; the high-PLL alpha is
`+0.11409501200` and `+0.11409501189`, respectively.  The corresponding
1e-10-to-1e-12 differences are 5.34×10⁻¹¹ and 1.09×10⁻¹⁰ s⁻¹.  The closest
7/8 pair agrees at similar or smaller differences.

The 1e-8 setting leaves equilibrium residuals around 3.52×10⁻⁹ and can
collapse a well-damped 7/8 alpha to approximately zero.  For V8 nominal it
also changes the sign from the converged `-0.00144` to approximately
`+1.38×10⁻⁸`.  This is a solver-accuracy failure mode, not evidence that the
converged V8 result is false.

For the ten key case/scenario combinations:

* `eigen`, `eigvals`, and complex-Float64 recomputation agree to at most
  2.11×10⁻¹² s⁻¹ in alpha;
* the descriptor generalized-eigenvalue cross-check differs by at most
  1.40×10⁻¹⁰ s⁻¹, with 89 finite non-gauge values for V8 and 92 for the
  7/8 cases;
* Float32 is not adequate here: its alpha error reaches 0.1441 s⁻¹ for a
  well-damped 7/8 case;
* BigFloat matrix eigensolution is unavailable in the installed Julia
  LAPACK interface and is explicitly marked unavailable in the raw table.

The reduced-Jacobian condition numbers in the key cases range from about
3.1×10¹⁸ to 8.4×10¹⁹, with smallest singular values as low as 6.1×10⁻¹⁴.
This ill-conditioning is present in 7/8 as well as V8, so it is a model-wide
numerical warning rather than a V8-only pathology.  It does require tight
equilibrium tolerances and argues against interpreting low-precision runs.

## Finite perturbations around V8

All 32 finite perturbation cases completed equilibrium and spectrum evaluation.
At nominal conditions, +/-0.5% PLL changes alpha from the baseline
`-0.0014448` to `-0.0032244` and `+0.0003819`; +/-0.5% load-P changes it to
`+0.0101162` and `-0.0103302`.  A +/-0.5% IBR-P change produces
`-0.0678135` and `+0.0788046`.  At the high-PLL corner V8 is already unstable;
the -0.5% IBR-P perturbation moves it to `-0.0141389` (stable but still below
the -0.05 engineering target), while +0.5% moves it to `+0.2548889`.

These are local finite-stress results, not a structured radius.  They reinforce
that V8 is close to the true-stability boundary nominally and sensitive to
operating point, but they do not identify a global perturbation boundary.

## Updated scientific judgment

* **Penetration versus composition:** the extreme 8/8 transition is primarily
  a cardinality/penetration endpoint in this case study; composition produces
  meaningful within-cardinality modulation and prevents a one-variable MW
  explanation.
* **Numerical status:** V8 is not dismissed as a numerical artifact.  Tight
  solver settings, independent Float64 eigensolvers, and the descriptor check
  reproduce the result.  The model is nevertheless ill-conditioned, so
  loose tolerances and Float32 are unacceptable.
* **Hidden dynamic margin:** this short audit did not run structured-radius
  searches or matched-alpha radius pairs.  It therefore neither confirms nor
  refutes the preregistered 3× hidden-fragility hypothesis.
* **Planning claim:** no new intervention or fresh holdout was run here.  The
  previous co-design and holdout conclusions remain bounded by their earlier
  report.

## Recommended next step

Keep the headline narrow: **a minimal 8/8 robustness blocker with
composition-sensitive margins and a documented conditioning hazard**.  Before
redesigning the IAS poster around hidden structured radius or universal
co-design value, run the explicitly deferred matched-alpha radius experiment
with the tight solver protocol and a predeclared cross-cardinality matching
rule.  If that radius result is absent, do not claim hidden fragility.
