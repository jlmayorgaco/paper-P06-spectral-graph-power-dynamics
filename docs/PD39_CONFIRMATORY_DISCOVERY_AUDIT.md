# PD39 confirmatory discovery audit

**Status:** completed before any confirmatory numerical experiment  
**Discovery commit audited:** `61554336`  
**Confirmatory branch:** `research/pd39-robust-transition-confirmatory`  
**Source worktree:** `C:\tmp\pd39-confirmatory`

This audit treats the first campaign as frozen discovery data. No model
evaluation was rerun to produce the statements below. The audit reads the
persisted CSV/TOML/Markdown outputs and the discovery source code.

## 1. What are the eight replacement candidates?

The candidate SG buses are:

`30, 32, 33, 34, 35, 36, 37, 38`.

Bus 31 is the slack bus and is excluded. A portfolio is a subset of these
eight buses, so the exhaustive campaign contains `2^8 = 256` portfolios.

## 2. What exactly are the nine scenarios?

The frozen discovery scenario set is the nominal case plus all eight corners
of three independent factors:

| Scenario | PLL factor | filter reactance factor Xf | current-control factor CC |
|---|---:|---:|---:|
| `nominal` | 1.0 | 1.0 | 1.0 |
| `pll0.8_xf0.8_cc0.8` | 0.8 | 0.8 | 0.8 |
| `pll0.8_xf0.8_cc1.2` | 0.8 | 0.8 | 1.2 |
| `pll0.8_xf1.2_cc0.8` | 0.8 | 1.2 | 0.8 |
| `pll0.8_xf1.2_cc1.2` | 0.8 | 1.2 | 1.2 |
| `pll1.2_xf0.8_cc0.8` | 1.2 | 0.8 | 0.8 |
| `pll1.2_xf0.8_cc1.2` | 1.2 | 0.8 | 1.2 |
| `pll1.2_xf1.2_cc0.8` | 1.2 | 1.2 | 0.8 |
| `pll1.2_xf1.2_cc1.2` | 1.2 | 1.2 | 1.2 |

The PLL gain parameters are regenerated consistently from the PLL bandwidth
factor. The discovery model is the stock PowerDynamics IEEE-39 network with
`PowerDynamics.Library.ComposableInverter.SimpleGFLDC` at replaced buses.

## 3. What exact mathematical quantity was called “robust margin”?

For a qualified equilibrium and its non-gauge eigenvalues, the discovery
code defines

```text
m_dyn(S,c) = - max Re(lambda_nontrivial(S,c)).
```

The portfolio-level robust margin is the minimum over the nine scenarios:

```text
m_9(S) = min_c m_dyn(S,c)
       = min_c [-max Re(lambda_nontrivial(S,c))].
```

Eigenvalues with absolute value at or below `1e-8` were treated as numerical
/ gauge modes and excluded from the engineering margin. The persisted
portfolio table contains `dynamic_margin = m_dyn(S,c)` and
`max_real = max Re(lambda_nontrivial(S,c))`.

## 4. Was the criterion `alpha <= -0.05 s^-1` or something different?

It is mathematically equivalent to that criterion, subject to the discovery
implementation's stored precision:

```text
dynamic_margin >= 0.05 s^-1
<=> -alpha >= 0.05 s^-1
<=> alpha <= -0.05 s^-1,
```

where `alpha = max Re(lambda_nontrivial)`. In code, robust feasibility was
implemented as `stable && dynamic_margin >= 0.05`; linear stability itself was
checked separately with the `1e-8` tolerance.

## 5. List all nine cases that failed the 0.05 target

The nine failures are exactly the nine V8 rows:

| Portfolio | Scenario | alpha [s^-1] | `m_dyn` [s^-1] | linearly stable? | robust target? |
|---|---|---:|---:|---|---|
| V8 | `nominal` | -0.001444778022811 | 0.001444778022811 | yes | no |
| V8 | `pll0.8_xf0.8_cc0.8` | -0.046447485569402 | 0.046447485569402 | yes | no |
| V8 | `pll0.8_xf0.8_cc1.2` | -0.044674976247513 | 0.044674976247513 | yes | no |
| V8 | `pll0.8_xf1.2_cc0.8` | -0.047001262373605 | 0.047001262373605 | yes | no |
| V8 | `pll0.8_xf1.2_cc1.2` | -0.045400154247959 | 0.045400154247959 | yes | no |
| V8 | `pll1.2_xf0.8_cc0.8` | 0.111414406182939 | -0.111414406182939 | **no** | no |
| V8 | `pll1.2_xf0.8_cc1.2` | 0.114095011886356 | -0.114095011886356 | **no** | no |
| V8 | `pll1.2_xf1.2_cc0.8` | 0.110443350865468 | -0.110443350865468 | **no** | no |
| V8 | `pll1.2_xf1.2_cc1.2` | 0.112906894696246 | -0.112906894696246 | **no** | no |

The worst V8 discovery robust margin is `-0.114095011886356 s^-1`.

## 6. Are they exactly V8 × nine scenarios?

Yes. The portfolio table contains `2304 = 256 × 9` rows and exactly nine
rows have `robust_feasible_for_scenario = false`. All nine have portfolio

```text
30;32;33;34;35;36;37;38
```

and all nine have `equilibrium_status = ok`. Thus these are not equilibrium
failures.

## 7. Do all 255 proper portfolios pass in all nine scenarios?

Yes. The persisted counts are:

```text
255 proper portfolios × 9 scenarios = 2295 robust-feasible rows.
```

The cardinality summary has robust-feasible counts `1, 8, 28, 56, 70, 56,
28, 8, 0` for cardinalities `0` through `8`, respectively. Therefore all
proper subsets, including every 7-of-8 predecessor, pass all nine scenarios.

## 8. Is V8 actually unstable (`alpha >= 0`) in any scenario?

Yes. V8 is linearly stable in the nominal case and in all four low-PLL
corners, but it is truly unstable in all four high-PLL corners:

```text
pll1.2_xf0.8_cc0.8
pll1.2_xf0.8_cc1.2
pll1.2_xf1.2_cc0.8
pll1.2_xf1.2_cc1.2
```

The positive rightmost real parts are listed in the table above.

## 9. Or is it nominally stable but insufficiently damped?

Both regimes occur. V8 is nominally stable but insufficiently damped in the
nominal condition and the four low-PLL corners. It is genuinely unstable in
the four high-PLL corners. Calling the whole V8 result merely “low margin”
would therefore be inaccurate.

For the nominal condition alone, V8 is not a true stability blocker; it is a
robustness blocker because its nominal margin is only
`0.001444778022811 s^-1`. For the nine-condition robust requirement, the
full V8 portfolio is both a true-stability blocker and a 0.05-robustness
blocker, after the confirmatory definitions below are frozen.

## 10. What are the nine exact alpha values for V8?

They are the `alpha` column in the table in Section 5. They are the exact
persisted `max_real` values printed to 15 decimal places; the source CSV
retains its full floating-point text representation.

## 11. Frequency, damping ratio, rightmost eigenvalue, critical mode family,
and dominant state/device participation for each V8 case

The discovery artifacts do **not** persist the raw complex eigenvalue vector,
eigenvectors, mode labels, frequencies, damping ratios, or participation
factors. They persist only the scalar rightmost real part (`max_real`) and
its negation (`dynamic_margin`), plus equilibrium and boolean status fields.

Consequently, the following quantities are not identifiable from discovery
data and are deliberately not reconstructed by inference:

| Quantity | Discovery status |
|---|---|
| frequency `abs(Im(lambda))/(2*pi)` | not archived |
| damping ratio | not archived |
| full rightmost complex eigenvalue | not archived; only `Re(lambda)` was archived |
| critical mode family | not archived |
| dominant state/device participation | not archived |

The confirmatory modal rerun must calculate and archive these quantities
before making a physical-mode claim. Until then, the only defensible
discovery statement is about the sign and scalar value of the rightmost
non-gauge real part.

## 12. Eight 7-of-8 predecessors

The following are the exact worst-case values across the nine discovery
scenarios. `alpha_worst = -m_9` because every predecessor is qualified and
stable in all nine scenarios.

| Missing SG bus | 7-of-8 portfolio | worst scenario | worst-case alpha [s^-1] | worst-case `m_9` [s^-1] |
|---:|---|---|---:|---:|
| 30 | `32;33;34;35;36;37;38` | `pll0.8_xf0.8_cc0.8` | -0.126228983466180 | 0.126228983466180 |
| 32 | `30;33;34;35;36;37;38` | `pll0.8_xf0.8_cc1.2` | -0.138022315452549 | 0.138022315452549 |
| 33 | `30;32;34;35;36;37;38` | `pll1.2_xf0.8_cc1.2` | -0.138295566944475 | 0.138295566944475 |
| 34 | `30;32;33;35;36;37;38` | `pll0.8_xf0.8_cc1.2` | -0.085834701546681 | 0.085834701546681 |
| 35 | `30;32;33;34;36;37;38` | `pll0.8_xf0.8_cc0.8` | -0.136356176223702 | 0.136356176223702 |
| 36 | `30;32;33;34;35;37;38` | `pll0.8_xf1.2_cc1.2` | -0.138057659072850 | 0.138057659072850 |
| 37 | `30;32;33;34;35;36;38` | `pll1.2_xf0.8_cc1.2` | -0.138325419498293 | 0.138325419498293 |
| 38 | `30;32;33;34;35;36;37` | `pll0.8_xf1.2_cc0.8` | -0.138019716389902 | 0.138019716389902 |

The required critical frequency, mode family, and participation for these
predecessors are likewise absent from the discovery scalar table and are
deferred to the preregistered confirmatory modal analysis.

## 13. Classification of the eleven outage equilibrium failures

The failed outage rows are:

| Branch | Endpoints |
|---:|---|
| 5 | 2–30 |
| 14 | 6–31 |
| 20 | 10–32 |
| 27 | 16–19 |
| 32 | 19–20 |
| 33 | 19–33 |
| 34 | 20–34 |
| 37 | 22–35 |
| 39 | 23–36 |
| 41 | 25–37 |
| 46 | 29–38 |

They are all stored only as `equilibrium_status = equilibrium_failed`, with
`dynamic_margin`, `max_real`, and `margin_loss` equal to `NaN`. The discovery
implementation's branch-removal wrapper catches the exception without
persisting its type or message. The saved table therefore cannot distinguish
islanding, graph disconnection, power-flow divergence, dynamic initialization
failure, numerical solver failure, or voltage infeasibility.

The correct audit classification is consequently:

> **unclassified / not estimable from discovery provenance** — not “unstable”
> and not a numerical category assigned after the fact.

The confirmatory campaign will rerun these eleven cases with typed diagnostic
stages (graph connectivity, power-flow convergence, finite PF state,
dynamic initialization, fixed-point residual, and spectrum) and will preserve
the original failure status alongside the new diagnostic provenance. Until
that rerun, assigning one of the requested physical causes would be
unsupported.

## Blocker definitions frozen for the confirmatory extension

For the nine-condition requirement, define:

```text
H_0.05 = minimal portfolios S such that
         every proper subset satisfies m_dyn(S',c) >= 0.05 for all c,
         while S violates that robust requirement.

H_0 = minimal portfolios S such that
      every proper subset is linearly stable for all c,
      while S has alpha(S,c) >= 0 for at least one c.
```

The discovery audit establishes the following conditional structure, pending
the confirmatory branch's frozen table implementation:

```text
H_0.05 = {V8}
H_0    = {V8}
```

These two statements are distinct. Under the nominal scenario alone, V8 is
stable and therefore is not a nominal true-stability blocker; its nominal
failure is only against the `0.05 s^-1` engineering margin.

## Audit verdict

The discovery campaign supports an exact finite-benchmark `255 + 1`
structure for the nine-condition robust screen. It also contains a true
stability failure for V8 in the four high-PLL corners. It does **not** yet
establish the physical mode mechanism, a structured dynamic radius, a
causal weak node/link, a 100% transition design, or a nonlinear TDS result.
Those are confirmatory questions and must be preregistered before new model
evaluations.

