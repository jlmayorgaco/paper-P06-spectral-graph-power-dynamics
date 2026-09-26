# E39 — baseline audit: the earlier AUC ≈ 0.14 reproduces, and sharpens

Every predictor is scored on its **raw value** and both the AUC and its
complement are reported, so no hidden sign convention can flatter or damn
anything. Sizes 4, 5 and 6 are audited separately because the base rate moves
from 0.08 to 0.93 across them and a pooled number would be meaningless.

Permutation sanity check: the label-permutation null mean over all 24 rows is
**0.50002**, as it must be.

## The lower-order reconstructions

| size | predictor | AUC (raw) | AUC (sign-flipped) | permutation p |
|---|---|---|---|---|
| 4 | additive | **0.144** | 0.856 | < 0.0002 |
| 4 | pairwise | **0.153** | 0.847 | < 0.0002 |
| 5 | additive | 0.146 | 0.854 | < 0.0002 |
| 5 | pairwise | 0.147 | 0.853 | < 0.0002 |
| 6 | additive | **0.090** | 0.910 | 0.0002 |
| 6 | pairwise | 0.105 | 0.895 | 0.0006 |

The earlier 0.14 / 0.15 headline **reproduces exactly** at size 4.

**Required wording.** These are **inversely ordered on this benchmark**. They are
*not* "universally anti-predictive": the complement is 0.85–0.91, which means the
reconstruction carries real ordering information with the wrong sign here, and
nothing in this experiment licenses a claim about other systems.

**The stronger and cleaner statement is about the threshold, not the ranking.**
Under the declared rule — predict unstable iff the reconstructed `α > 0` — the
lower-order reconstructions **never fire at all**:

| size | TP | FP | TN | FN | recall | balanced accuracy |
|---|---|---|---|---|---|---|
| 4 | 0 | 0 | 116 | 10 | 0.000 | 0.500 |
| 5 | 0 | 0 | 75 | 51 | 0.000 | 0.500 |
| 6 | 0 | 0 | 6 | 78 | 0.000 | 0.500 |

Every reconstructed value is at or below zero. A model that keeps only
interactions up to order 3 does not merely rank these portfolios badly — **it
declares every one of them stable**, including the 139 that are not. That is the
Track-A claim stated as a classifier property, and it is the form that belongs in
a paper.

## The conventional predictors, fairly reported

| predictor | AUC in the conventional direction | verdict |
|---|---|---|
| removed inertia fraction | 0.87 / 0.78 / 0.85 (sizes 4/5/6) | **the best conventional ranker** |
| replaced MW | 0.81 / 0.77 / 0.86 | **good** |
| max MIIF | 0.53 / 0.55 / 0.75 | weak to moderate |
| min SCR (low = risky) | 0.60 / 0.65 / 0.68 | weak |
| mean SCR (low = risky) | 0.59 / 0.62 / 0.66 | weak |
| gSCR (low = risky) | 0.64 / 0.63 / 0.72 | weak |

The short-circuit family reproduces E13's published numbers once the sign
convention is made explicit: E13's `min_scr` AUC of 0.605 at size 4 is exactly
the complement of the raw-value 0.395 reported here.

**This corrects an overstatement.** "Conventional baselines are anti-predictive"
is not supportable. Removed inertia and replaced megawatts rank instability
reasonably well, AUC 0.77–0.87. What they cannot do is say *which particular*
portfolio of a given size fails — precision at 10 is 0.30–1.00 and is driven by
the base rate — and they carry no mechanism. The defensible claim is that the
short-circuit family is weak here and that no conventional measure identifies the
irreducible fourth-order structure.

## Files

`E39_baseline_audit.csv`, `E39_ROC_PR.png`, `E39_figure_source.csv`,
`manifest.json`.
