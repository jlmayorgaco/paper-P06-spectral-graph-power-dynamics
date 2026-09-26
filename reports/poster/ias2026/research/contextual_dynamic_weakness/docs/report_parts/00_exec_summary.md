# Contextual Dynamic Weakness in Inverter-Rich Power Networks

**CDW campaign final report.** Branch `research/contextual-dynamic-weakness`,
created from tag `TX4_FINAL_MANUSCRIPT_FREEZE` (commit `69f200df`). Not
pushed. TX4 was read only and never modified.

## 1. Executive summary

The campaign tested the meta-hypothesis

> static weakness ≠ dynamic weakness ≠ contextual weakness ≠ intervention
> leverage

against the frozen TX4 IEEE-39 model, using only genuinely new preregistered
experiments (never rerunning or reinterpreting any TX4 result).

**It does not hold as an equality anywhere, and it does not collapse to a
single ranking anywhere either — every distinction the hypothesis draws is
empirically real on this benchmark, at material effect sizes.**

**GOLD gates:** **A and B PASS**; **C, D, E and F do not** (E/F are BLOCKED,
not failed, by a missing-material gate). Fifteen of twenty-four testable
sub-claims are SUPPORTED, eight are NOT_SUPPORTED (seven of these are
informative negatives with an identified mechanism, not inconclusive
failures), one is deferred by design.

**Strongest results.**
- **Contextual sign reversal is real, recurrent, and survives the strongest
  available scrutiny** (mode tracking, an oracle-fitted node ranking, and four
  independent uncertainty envelopes) — **GOLD-A**.
- **Total, re-equilibrated dynamic line sensitivity beats every static
  baseline by a wide, cross-validated margin** (0.987 vs 0.50 median Spearman
  against held-out finite interventions) — **GOLD-B**, the single strongest
  result of the campaign.
- **Topology alone, with no control retuning, both removes and creates
  transverse incompatibilities**, and no static topology score anticipates
  which.
- **Weak corridors exist and are structurally, not spectrally, defined** —
  transformer groups beat every graph-spectral partition tested.
- **A genuinely reduced model is certifiable and accurate, but buys no
  end-to-end speed** on this benchmark, because the equilibrium solve, not
  the eigenanalysis, is the bottleneck — a precise, actionable negative
  result (**GOLD-C: FAIL**).
- **Plan-level multi-constraint design never had to demonstrate an advantage**
  at the scale it could be cleanly tested: single-boundary tuning already
  sufficed in every case (**GOLD-D: FAIL**, informatively).
- **The Africano/PV benchmark is blocked** for lack of source material
  (**GOLD-E/F: BLOCKED**), exactly as the preregistration required rather than
  substituting a different feeder.

**No weakness index was constructed anywhere**, consistent with the
preregistration and reinforced by the campaign's own evidence (E10: an
averaged attribution is dominated by rare fast-mode contexts and hides the
real, modest electromechanical reversals).

**Determinism.** Every rerun deterministic case (E1 census at three policies,
E4 link sensitivities at P4, E9 at target T1, all families/methods) was
byte-identical on its summary fields, excluding wall-clock.
