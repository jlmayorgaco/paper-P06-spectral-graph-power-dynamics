# F1B — excitation provenance audit: **the phenomenon survives, the "order 4" label does not**

**Verdict.** The higher-order phenomenon is **not** an artifact of the inherited
excitation point. It occupies a finite two-dimensional region of the **base-stable,
model-plausible** excitation parameter region, and it occurs at points where **both** excitation parameters lie
inside the documented ranges of the model family. The IEEE-39 flagship is
therefore **not** downgraded on the ad-hoc-parameter criterion.

But the audit forces a different and larger correction: **the composability order
is not a property of the benchmark. It is a function of the excitation dynamics**,
taking every value in {1, 2, 3, 4, ∞} across the base-stable swept region. "Order 4" is
one cell of a map, not a fact about IEEE-39.

## What the inherited parameters are, and are not

The internal AVR is

    efd' = (K (vref + vs - |V|) - efd) / T

with `K = KA = 10.1` and `T = TE` taken **per machine** from the case, where TE
ranges over 0.25–0.50 s. These were obtained by pairing the IEEEX1 **amplifier**
gain with the IEEEX1 **exciter** time constant — two different blocks — so the
model is a substitution, not a reduction, and was not entitled to be called
representative before this audit.

Provenance reference: the **ANDES SEXS model card**, the reference implementation
of this first-order family. With the lead-lag collapsed (`TA/TB = 1`) SEXS is
exactly `K/(1 + sTE)`. Its documented metadata:

| parameter | default | documented range | inherited value |
|---|---|---|---|
| `K` | 20 | **[20, 100]** | **10.1 — below the range** |
| `TE` | 1 | **[0, 0.5]** | 0.25–0.50 — **inside** |
| design note in the source | | `5 ≤ K·TA/TB ≤ 15` | 10.1 — **inside** |

So the inherited gain is below the stated `K` range while satisfying the design
guideline, and the time constants are inside their range. **No single point
satisfies all three documented criteria at once**, and that is a property of the
guidance, not of this study: the vrange minimum of 20 and the design maximum of
15 do not overlap when `TA/TB = 1`.

## Two sweeps, because uniformity is itself a model change

**Flattening the fleet destroys the phenomenon.** At `K = 10.1` with a *uniform*
`T = 0.25 s` on every machine the flagship is **stable at −0.197, κ = ∞**. With
the case's own heterogeneous `TE` at the same gain the flagship is **unstable at
+0.145, κ = 4**. Real fleets have diverse exciters, so the multiplicative sweep —
which scales every machine's own value and preserves the diversity — is the
physically meaningful instrument, and the uniform grid is reported beside it only
to compare against documented ranges.

Branch-independent order used throughout, the definition F2 formalises:

    Gamma = right half plane intersected with the frozen 0.3-1.5 Hz band
    N(S)  = modes of the replaced system inside Gamma
    kappa = min { |S| : N(S) != N(empty) },   infinity if no subset qualifies

No modal labels appear anywhere in it.

### Multiplicative sweep, PSS on — 38 admissible of 72

| `K/K_case` \ `T/T_case` | 0.25 | 0.50 | 0.75 | 1.00 | 1.25 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|
| 0.50 | **4** | ∞ | ∞ | ∞ | ∞ | **4** | ∞ | ∞ |
| 0.75 | — | 2 | ∞ | ∞ | **4** | **4** | 3 | ∞ |
| **1.00** | — | — | 3 | **4 ← nominal** | **4** | **4** | 3 | 3 |
| 1.25 | — | — | — | 3 | 3 | 3 | 3 | 2 |
| 1.50 | — | — | — | 1 | 3 | 3 | 3 | 2 |
| 2.00 | — | — | — | — | — | 2 | 3 | 2 |
| 2.50 | — | — | — | — | — | — | 2 | 2 |
| 3.00 | — | — | — | — | — | — | — | 2 |
| 4.00 | — | — | — | — | — | — | — | 1 |

(— = base case unstable, rejected per item 6.)

`κ = 4` at **7 of 38** base-stable points, forming a region around the nominal
rather than a knife edge. Distribution: `∞`:9, 1:2, 2:8, 3:12, **4:7**.
PSS off is very similar: `κ = 4` at 6 of 38.

"Base-stable" is the honest qualifier throughout: see the wording constraint at
the end of this document.

### Uniform sweep — the documented-range test

Of the 25 base-stable points with `K` inside the documented [20, 100],
**6 give κ = 4**. The cleanest is

> **K = 20, T = 0.40 s** — gain inside [20, 100], time constant inside [0, 0.5],
> base stable at −0.049 (PSS on) / −0.078 (PSS off), flagship unstable at +0.072 /
> +0.081, **κ = 4** in both stabilizer settings.

This is a fully documented-admissible `κ = 4` point, obtained without tuning
anything toward instability. It is the answer to item 7: the phenomenon occupies
a finite region of model-plausible parameters, not only the inherited point.

## The network is fragile to excitation gain

**34 of 72** points in the multiplicative sweep and **42 of 84** in the uniform
sweep were rejected because the **base case** is unstable — no replacement
involved. The largest effective gain that leaves the base stable is about 40.
That is a property of this network with a first-order exciter and it sharply
limits the admissible region. It is also why the distributed IEEEX1 data, with
`KE = −0.05`, gives six unstable real modes before anything is replaced.

## An unplanned finding: the closure margin orders the composability order

Median flagship closure margin over the base-stable multiplicative sweep:

| κ | points | median `m_cl` |
|---|---|---|
| 1 | 2 | 0.359 |
| 2 | 8 | 0.340 |
| 3 | 12 | 0.237 |
| **4** | **7** | **0.073** |
| ∞ | 9 | 0.166 |

The closure margin at the full portfolio falls monotonically as the composability
order rises: the harder the failure is for a lower-order model to see, the closer
the interaction operator sits to its singularity. This is **descriptive** — it was
not preregistered and it is one benchmark — but it is a coherent link between the
combinatorial and the operator-theoretic sides of the work, and it is the natural
thing for F2/F5 to try to prove.

## Consequences for the claims

**Survives.** The higher-order phenomenon is real over a finite, base-stable
region that includes points inside the documented parameter ranges. The IEEE-39 flagship is not downgraded for
ad-hoc parameters.

**Must be corrected.** "IEEE-39 exhibits a fourth-order non-composability" is not
supportable as stated. The supportable claim is:

> On IEEE-39 with a first-order excitation model, the spectral composability
> order of the 30+33+35+37 portfolio takes values in {1, 2, 3, 4, ∞} depending on
> the excitation gain and time constants, and equals 4 over a region that
> contains the case's own excitation data and points inside the documented
> parameter ranges of the model family.

**Strengthened.** F1 showed composability depends on excitation *model*; F1B
shows it depends continuously on excitation *parameters*, and that the order can
be moved to 1, 2, 3, 4 or ∞ within the base-stable swept region. Together with E30's
reactive-policy result this is the substantive finding: **spectral composability
order is a policy-dependent property of network, dynamics and control together.**

**Still required on every claim.** The reactive dispatch (E30), the machine data
(E37), the excitation model (F1) and now the excitation parameters (F1B).

## Files

- `docs/F1B_EXCITER_PROVENANCE.md` — this document
- `results/F1B_exciter_parameter_map.csv` — 312 points, both sweeps, both stabilizer settings
- `figures/F1B_exciter_stability_map.png`, `figures/F1B_figure_source.csv`

## Wording constraint carried forward

The SEXS model card contains an internal tension: with the lead-lag collapsed
(`TA/TB = 1`) the stated gain range `K in [20, 100]` and the design note
`5 <= K*TA/TB <= 15` **cannot both be satisfied**. Until a literature audit
resolves which is authoritative for this family, the region swept here is
described as the **base-stable, model-plausible excitation region** and never as
"the entire physically admissible region". The ambiguity is a property of the
guidance, not of this study, and it is not hidden.
