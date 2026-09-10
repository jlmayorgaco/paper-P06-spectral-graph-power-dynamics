# E38 — is bus 30 individually ordinary but collectively enabling? **NOT SUPPORTED AFTER ADJUSTMENT**

The protocol's own rule decides this: *do not claim it if bus 30 loses
significance after controlling for megawatts and inertia*. It does. The claim is
therefore **not made**, and the evidence on both sides is recorded because it is
genuinely two-sided.

512 portfolios from the E12 census, 9 candidate buses, 42 minimal incompatible
cores (unstable with every proper subset stable), of which **12 are genuine
inter-area cores** and 153 portfolios fail in the PLL family.

## Bus 30 is individually ordinary — confirmed

| measure | bus 30 | rank among 9 |
|---|---|---|
| damping when replaced alone | −0.1390 | **5th** |
| participation in the base inter-area mode | 0.0265 | **7th** (near the lowest) |
| megawatts replaced alone | 1040 | 4th |
| mean minimum SCR when present | 1.630 | tied highest |
| share of all minimal cores | 0.524 | **6th** (bus 33 leads at 0.690) |

On every individual measure bus 30 is unremarkable, and on mode participation it
is nearly the *least* involved machine in the base inter-area mode.

## The unadjusted collective evidence is strong

| test | result |
|---|---|
| **share of the 12 genuine inter-area cores containing bus 30** | **12 of 12, 100 %** — the only bus at 1.00 |
| Fisher exact, order-4 portfolios, 9/56 unstable with bus 30 against 1/70 without | odds ratio **13.2**, **p = 0.0050** |
| matched pairs sharing three members, bus 30 swapped for another bus | **41 wins, 0 losses, 239 ties**, McNemar exact **p = 9.1e-13** |
| label permutation, 20000 draws | difference 0.146, **p = 0.0026** |

Bus 30 never loses a matched comparison. In 41 discordant pairs of order-4
portfolios that differ only in whether bus 30 or some other bus joins the same
three machines, the bus-30 portfolio is the unstable one every time.

## The adjusted test does not resolve it

Logistic regression on the 126 order-4 portfolios,
`unstable ~ bus30 + replaced_MW + removed_inertia + min_SCR + gSCR`:

| | value |
|---|---|
| bus-30 odds ratio, adjusted | **4.09** |
| 95 % interval | **[0.14, 124.2]** |
| p | **0.418** |

**Not significant, and the interval spans three orders of magnitude.** With 10
unstable cases and 5 predictors this model has almost no power; the restricted
inter-area-only fit did not converge at all (complete separation) and its numbers
are not reported. So the correct reading is *the adjustment cannot resolve
whether bus 30 matters beyond megawatts and inertia*, not *megawatts and inertia
explain bus 30 away*.

That distinction does not rescue the claim. The frozen rule was stated in advance
and it is not met.

## Status

- **Claimable:** bus 30 appears in every one of the 12 genuine inter-area minimal
  cores while ranking 5th, 6th and 7th on the individual measures. That is a
  *descriptive* fact about this census and may be stated as one.
- **Not claimable:** that bus 30 is a collective structural enabler distinct from
  its megawatts and inertia. The multivariable test does not support it.
- **What would settle it:** a larger unstable sample — the held-out Monte Carlo
  E35 produces 240 unstable operating points, and repeating this audit over the
  per-sample lattices there would give the adjustment real power. Recorded in
  `TPWRS_REMAINING_GAPS.md`.

A methodological correction was made during the run: the matched-pair test was
first computed with a 2×2 Fisher table, which is the wrong test for paired data.
It is now McNemar's exact binomial on the discordant pairs. The conclusion is
unchanged.

## Files

`E38_bus30_collective_role.csv`, `E38_bus30_logistic.csv`, `manifest.json`.
