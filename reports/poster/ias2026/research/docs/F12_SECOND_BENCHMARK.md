# F12 — Second documented benchmark: Kundur two-area

**Verdict under the preregistered rule: REPRODUCED on both slices** — physical
reactive policy changes the minimal incompatibility structure of the same
replacement fleet on an independently parameterized grid, through imaginary-axis
crossings detected by the port closure. **But the reproduction carries a strong
qualifier that was not anticipated in the preregistration**: on Kundur, every
pair of replacements carries an aperiodic (zero-frequency) right-half-plane mode
that the 0.3–1.5 Hz window excludes by construction, so the band-limited
hypergraph there is not the operative safety object and the fleet is never fully
composable.

## Protocol

Preregistered and committed before any F12 composability number
(`configs/kundur/F12_preregistration.yaml`, sha256 `a545fc7e…`, commit
`34a64bdc`). IEEE 68-bus/NETS-NYPS was preferred but no documented dynamic
dataset is available offline here; the Kundur case bundled with ANDES 2.0.0
(`kundur_sexs.xlsx`, sha256 `7ad5a575…`) was used.

- Documented models: GENROU (Kundur Example 12.6 data), SEXS (`K = 20`,
  `TE = 0.83`, `TA/TB = 0.4`, `TB = 5`), TGOV1, no PSS. Declared reductions, the
  same as IEEE-39: two-axis machines, no governor, constant-power loads; SEXS
  represented exactly (a lead-lag exciter state was added to the model for this,
  default-off for IEEE-39). No parameter changed.
- Eligibility checks: Ybus against ANDES `1.8e-15`; power flow `2.3e-7`; base
  stable (`-0.122`); inter-area mode `-0.145 +- j 2pi 0.657 Hz` (ours) versus
  `-0.154 +- j 2pi 0.650 Hz` (ANDES full model).
- Candidates `{2, 3, 4}` (every non-slack machine, IEEE-39 rule); the F7 converter
  and leaky Q/V policy; band 0.3–1.5 Hz; maps K12A (`g x k`) and K12B (`g x t`),
  61 x 61 nodes each; edge bisection, Safeguard B, port closure.

## Results

| | K12A (`g x k`) | K12B (`g x t`) |
|---|---|---|
| base-unstable nodes | 0 | 1 525 of 3 721 (`t > 1.4`) |
| `kappa` values | 1, 2, inf | 1, 2, inf |
| distinct `H` (band) | 9 | 7 |
| pure-policy lines on which `H` changes | **61 of 61** | **36 of 36** |
| located crossings | 554 | 327 |
| imaginary-axis / lower-band-edge | 442 / 112 | 245 / 82 |
| imaginary-axis crossings port-visible | 442 / 442 (closure at most `1.1e-6`) | 245 / 245 (at most `1.2e-6`) |
| pure-policy imaginary-axis `H` changes | 327 | 182 |
| random off-grid checks mismatching the nearest node | 9 / 200 | 2 / 200 |

A typical pure-policy line (K12A, `k = 0.75`), as `g` rises from 0:

    EMPTY -> {2,4} -> {2,3 | 2,4} -> {2} -> {2,4} -> EMPTY      (kappa inf -> 2 -> 2 -> 1 -> 2 -> inf)

so, as on IEEE-39, voltage regulation at intermediate gain creates
incompatibility that neither fixed-Q nor full regulation shows. Crossing
frequencies: witness `{2}` 0.31–0.73 Hz (median 0.63, the inter-area range);
`{2,3}` 0.75–0.82 Hz; `{2,4}` 0.805 Hz; `{3}`, `{4}` 0.30–0.43 Hz, near the lower
band edge.

## The qualifier: outside-band instability

| | points with `H = EMPTY` | of which some subset has an RHP mode outside the band |
|---|---|---|
| IEEE-39 F7A / F7B / F7C | 12 041 / 11 208 / 19 187 | **0 / 0 / 0** |
| Kundur K12A / K12B | 556 / 374 | **556 / 374** (all) |

On Kundur every pair `{2,3}`, `{2,4}` (and `{3,4}` at 93 % of those points)
diverges aperiodically — a real eigenvalue between `+20` and `+1 250 s^-1`, the
grid-following loss of synchronism in a grid too weak to host two of four units
as current sources (the same mechanism as IEEE-39 Track B, N6c). Singles at fixed
Q additionally carry 0.2 Hz oscillations just below the band.

**Post-hoc sensitivity (declared as such, not preregistered):** recomputing `H`
with the protected region equal to the whole right half plane
(`results/F12/F12_full_rhp_sensitivity.json`):

| | full-RHP `H` equal to band `H` | distinct full-RHP `H` | policy lines on which full-RHP `H` changes | points with full-RHP `H = EMPTY` |
|---|---|---|---|---|
| IEEE-39 F7A / F7B / F7C | 100 % / 96.5 % / 100 % | 30 / 36 / 16 | — | as in the band maps |
| Kundur K12A | 8 % | 4 | 61 of 61 | **0** |
| Kundur K12B | 0.3 % | 5 | 36 of 36 | **0** |

On IEEE-39 the band-limited hypergraph *is* the full-stability hypergraph. On
Kundur the full-RHP hypergraph is still policy-dependent (it moves between
`{2 | 3,4}`, `{2,3 | 2,4 | 3,4}`, `{2 | 3 | 4}` and `{2 | 4}`), but it is never
empty: some coalition of two replacements is always unsafe, whatever the policy.

## What F12 establishes

- **Framework-general:** the incompatibility hypergraph is policy-dependent on a
  second, independently parameterized grid, under both the preregistered and the
  full-RHP protected region; its changes are imaginary-axis crossings where the
  port closure reaches `-1`; voltage regulation is non-monotone there too.
- **Benchmark-specific:** `kappa = 4`, the flagship coalition, the tongue, and a
  safe region reachable by policy are IEEE-39 facts. Kundur has `kappa <= 2`
  (three candidates) and no policy makes its fleet fully composable.
- **Methodological lesson:** a band-limited protected region is only as good as
  its band. Where aperiodic instability is possible — weak grids with
  grid-following converters — `Gamma` must include the real axis, or the band
  result must be reported beside the full-RHP one, as here.

Files: `results/F12/` (nodes, events, spot-check, outside-band and full-RHP
analyses), `figures/F12_kundur_maps.png`.
