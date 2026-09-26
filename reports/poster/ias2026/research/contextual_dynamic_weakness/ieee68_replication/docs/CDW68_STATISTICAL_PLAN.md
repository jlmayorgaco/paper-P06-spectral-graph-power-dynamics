# CDW68 statistical plan

Companion to `docs/CDW68_PREREG_V1.md`. It is frozen with the preregistration.

## 1. Nature of the data

- Policies, draws and branch conditions are designed points (Latin hypercube; declared windows and envelopes). They are not samples from a physical population.
- Within a condition, the 64 portfolios are a complete census.
- Every fraction is a deterministic **coverage** of the tested conditions and is never called a probability.
- Effect sizes come first; p-values are secondary summaries under an exchangeability assumption stated for each test.

## 2. Reversal (R7, R13, R15)

- **Coverage:** k/n over eligible (base-stable) conditions.
- **Intervals:** percentile bootstrap over conditions (10^4 resamples, seed 20260920) and the exact Clopper–Pearson interval. Both are descriptive.
- **Gate:** it is decided on the point estimate. The intervals are reported next to it.
- **Magnitudes:** quartiles and maximum of min(|Δ(S1)|, |Δ(S2)|) over pairs, together with the counts of distinct (policy, unit) and distinct S2 marginals. Pairs overlap, so pair counts are not evidence weights.
- **Damping-ratio equivalent:** |Δσ| / ω of the tracked mode.

## 3. Fixed ranking (R9)

- **Per-policy measures:** p* and regret; medians and IQR over base-stable policies, with bootstrap intervals.
- **T1:** exact one-sided sign test of p* < 0.90 over eligible REAL Model A policies. Exchangeability of policies under the null "p* ≥ 0.90 is as likely as not" is assumed.
- **Transfer:** off-diagonal against diagonal regret, reported with the increment.

## 4. Branch ranking (R11, R12)

Per condition, against the finite truth at each γ:
- Spearman and Kendall;
- concordance;
- top-3 and top-5 precision;
- NDCG@5;
- sign accuracy on branches with |finite| ≥ τ.

Rules:
- **Gate statistics:** medians over the 12 holdout conditions.
- **Paired advantage:** ρ(D_tot) − ρ(comparator), per condition. Primary: the median and its percentile bootstrap interval over conditions (10^4, seed 20260921).
- **T2A / T2B:** one-sided sign-flip test of the median (advantage − 0.20) over conditions (exact enumeration, 2^12). Wilcoxon is reported as well.
- **Holm family** F68 = {T1, T2A, T2B}: raw and Holm-adjusted p-values are both reported.
- **R11:** sign agreement, median |D − FD| / max(|FD|, 1e-4), and Spearman over the 80 (branch × policy) pairs per model.

## 5. TDS (R14)

- Descriptive only.
- For each case: the linear Δα, the TDS Δσ̂ in each context, sign agreement, and the run label (RECOVERS / FAILS / OUTSIDE_MODEL_SCOPE / NUMERICAL_FAILURE).

## 6. Uncertainty (R15)

- **Phenomenon robustness:** the fraction of draws with at least one reversal (level D, and level T), per envelope.
- **Witness identity:** the fraction of draws whose strongest witness equals the nominal one.
- **Ranking:** Spearman quantiles over the ranking subset.
- No probability language.

## 7. Multiplicity and reporting

- Only F68 is Holm-corrected. Every other quantity is descriptive.
- Negative results get the same prominence as positive ones.
- Raw per-condition values are written to `results/`.
