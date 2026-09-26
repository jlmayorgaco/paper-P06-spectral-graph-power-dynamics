# CDW hardening — statistical plan (frozen with the preregistration)

## 1. What the data are

- **Policy points** are designed. Neither the old holdout nor the new
  HARDENING_H01–H24 are random samples of a physical population.
- **Envelope draws** are Latin-hypercube draws from declared bounds. They are not
  measured physical uncertainty.
- Every interval below is therefore a **resampling interval over the tested
  conditions**. It is a statement about the stability of a summary across this
  design, not a population confidence statement.

Required wording:
- "coverage across the tested holdout policies is X";
- "coverage fraction over the declared envelope draws is X";
- never "the probability of reversal is X".

## 2. Effect sizes first

Every result is reported as:
- the raw per-condition values, in the supplement CSV;
- the median and IQR (plus the mean where stated);
- a bootstrap interval.

A p-value is secondary. It is shown only with an explicit null and exchangeability
statement, as given in §4.

## 3. Bootstrap

- **Unit.** The resampling unit is the cluster:
  - GOLD-A: one policy = one cluster.
  - GOLD-B primary set: the 24 HARDENING policies are 24 clusters, and the 40
    fresh draws form 4 clusters by envelope (28 clusters). Within a resampled
    cluster, all of its conditions enter.
  - Corridors: as for GOLD-B.
- **Procedure.** B = 10 000 resamples, seed 20260932. The interval is the 2.5–97.5
  percentile (two-sided), and 5th percentile for one-sided lower bounds. The
  statistic is recomputed on the pooled conditions of each resample.
- **The H7 gate** uses the two-sided 95 % interval's lower end (2.5th percentile)
  for the median Δρ. This is conservative relative to a one-sided bound.
- **Old vs new shift.** An independent bootstrap of each set; the interval is for
  the difference of medians.

## 4. Inferential summaries (secondary) and their nulls

| test | statistic | null | exchangeability / interpretation |
|---|---|---|---|
| T1 (GOLD-A, new holdout) | number k of base-stable policies with p* < 0.90; exact one-sided binomial sign test | median p* ≥ 0.90 over the design | Conditional on the design, the policy-level outcomes are treated as exchangeable signs about 0.90. Descriptive of the design, not of a population |
| T2 (GOLD-B, H7) | median of Δρ − 0.20; one-sided cluster sign-flip permutation, 10⁴ flips, seed 20260933 | the cluster-mean differences are symmetric about 0.20 | Sign-symmetry of paired differences under exchangeable clusters; a secondary summary of the paired effect |
| T3 (mixing, H13 primary) | Spearman(μ(θ, ω_ref), n_rev); two-sided permutation of n_rev over policies, 10⁴, seed 20260933 | no monotone association across policies | The policies are exchangeable labels under the null |
| secondary (F_mix) | as T3, for the other predictors and targets | same | Holm within F_mix |
| descriptive (H4) | Mantel permutation of distance vs degradation | no association | Descriptive only |

## 5. Multiple testing

- **Confirmatory family F_conf = {T1, T2, T3}**, with Holm correction at family-wise
  α = 0.05. T3 additionally needs a Holm-adjusted p ≤ 0.01 and abs(ρ) ≥ 0.5 to
  "correlate", consistent with the old E12 bar.
- **Gates are effect-size thresholds** (§2 of the prereg). No gate is converted into
  or replaced by a p-value.
- **Secondary family F_mix.** Holm within the family. These results are labelled
  secondary.
- **Everything else is exploratory** and labelled so. Confirmatory language
  ("confirms", "holds on the new holdout") is used only for the prereg gates and
  for F_conf.

## 6. Correlations

- **Measures.** Spearman is primary; Kendall τ_b is reported with it. Ties use the
  scipy defaults.
- **Minimum sample.** A correlation over fewer than 6 finite pairs is reported as
  NaN and not interpreted.

## 7. Rank metrics (GOLD-B)

- **Orientation.** Truth t_e = −Δα_e, so larger means more stabilizing. Every
  predictor is oriented the same way (old E4 orientation).
- **Pairwise concordance** is the fraction of pairs with abs(t_e − t_f) ≥ τ_res that
  the predictor orders the same way. Predictor ties count ½.
- **Top-k precision** is abs(top-k(pred) ∩ top-k(truth))/k.
- **NDCG@5** is Σ_{r≤5} rel(π_r)/log₂(r+1) divided by the same sum for the ideal
  order, with rel_e = t_e − min_f t_f.
- **Sign accuracy** covers branches with abs(t_e) ≥ τ_mat, and compares the predictor
  sign with the truth sign. For static predictors, which are nonnegative, sign
  accuracy is not defined and is reported as NaN.

## 8. Handling non-finite and failed cases

- Failed solves are NaN and are counted. Metrics use the finite pairs.
- A condition with fewer than 20 finite branches is ineligible and is counted.
- No imputation.

## 9. Stratification (H26)

Every headline result is reported in three strata:
- FULL;
- EM-only (0.1–2.0 Hz);
- same-mode tracked, where meaningful.

A result "survives fast-mode removal" if its gate or effect direction holds in the
EM stratum.
