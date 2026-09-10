# F1 — ANDES model reconciliation: **HARMONIZED_AGREE**

**The two implementations were never disagreeing about physics. They were solving
different excitation models.** With an equation-equivalent AVR, the internal
simulator and ANDES 2.0.0 agree across the entire electromechanical band, mode
for mode, to **1.6e-6**, including the unstable flagship.

The G2 finding from the overnight run is **superseded**.

## The ladder, and where the sign flipped

Progressive harmonization, all with the governor removed and constant-power loads
on both sides. The observable is the worst mode in the frozen 0.3–1.5 Hz band.

| stage | what is harmonized | internal flagship | ANDES flagship | agree? |
|---|---|---|---|---|
| **R0/R1** | network, shunts, loads, equilibrium | — | — | **yes**, V to 1.7e-7 pu, generator S to 2.0e-6 pu |
| **R2** | + manual excitation, no stabilizer | **−0.25417 @ 0.9327 Hz** | **−0.25417 @ 0.9327 Hz** | **yes**, α to 1.0e-7 |
| **R3** | + first-order AVR, no stabilizer | **+0.32873 @ 0.5825 Hz, RHP 1** | **+0.32873 @ 0.5825 Hz, RHP 1** | **yes**, α to 1.5e-6 |

Band-restricted spectrum comparison, all twelve reconciliation cases:

| | value |
|---|---|
| worst band eigenvalue matching distance | **1.618e-06** |
| band mode **count** agrees | **12 of 12** (9,8,7,8,8,5 at R2; 11,10,10,10,10,7 at R3) |
| band right-half-plane count agrees | **12 of 12** |
| worst band α error | 1.5e-6 |
| worst band frequency error | 1.3e-7 Hz |

The *full* spectra do **not** match, and should not: ANDES's degraded GENROU
retains damper-winding modes near 1e4 rad/s that a fourth-order machine does not
have, so the whole-spectrum nearest-neighbour distance is 0.13–0.42. Every mode
inside the electromechanical band is matched to 1.6e-6. That is the invariant the
entire argument depends on.

## The subsystem that caused the disagreement

Service ablation on **both** sides, static-injection configuration:

| configuration | internal: base → flagship | ANDES: base → flagship |
|---|---|---|
| AVR + PSS | −0.126 → **+0.262 destabilizing** | −0.206 → −0.329 stabilizing |
| AVR, no PSS | −0.160 → **+0.329 destabilizing** | −0.205 → −0.328 stabilizing |
| field frozen, no PSS | −0.186 → **−0.260 stabilizing** | −0.216 → −0.338 stabilizing |

**The AVR is the first subsystem that changes the sign disagreement.** ANDES is
almost completely insensitive to its own AVR and PSS; the internal model is not.

The cause is a parameter pairing. The internal AVR is

    efd' = (K (vref + vs - |V|) - efd) / T,   K = case KA = 10.1,  T = case TE = 0.25

which takes the **amplifier gain** and the **exciter time constant** from two
different blocks of IEEEX1 and fuses them into one lag. IEEEX1 is an amplifier
`KA/(1+sTA)` with `TA = 0.06`, an exciter `1/(sTE + KE)` with `KE = −0.05`, and
rate feedback `KF1 s/(1+sTF1)`. The internal surrogate is therefore a
**substitution**, not a reduction — a faster, higher-gain excitation loop, which
is classically destabilizing for inter-area modes.

Once ANDES is given that same first-order AVR — SEXS with `TA/TB = 1`, which
collapses its lead-lag to unity and leaves `K/(1+sTE)` exactly — **ANDES produces
the unstable flagship too, at the same eigenvalue.**

## The consequence that matters more than the reconciliation

The flagship instability is **conditional on the excitation model substituted for
the benchmark's documented exciter**:

| excitation | flagship |
|---|---|
| first-order AVR, K = 10.1, T = 0.25 (harmonized) | **+0.329, unstable** — both tools |
| manual excitation (field frozen) | **−0.254, stable** — both tools |
| IEEEX1 as distributed | **−0.329, stable** — ANDES |

And the distributed IEEEX1 data cannot simply be adopted instead: with
`KE = −0.05` it gives the **base case six unstable real modes** and a spectral
abscissa of `+1.03` before any replacement. The IEEE-39 dynamic dataset as
distributed is not a usable small-signal stability benchmark, which is why a
substitution was made in the first place.

**So the honest framing is:** Track A is a property of the IEEE-39 *network* plus
a standard first-order excitation model. It is not a property of the IEEE-39
benchmark as distributed, and it is not robust to the choice of excitation model.
Every claim must carry that condition, exactly as it must carry the reactive
dispatch (E30) and the machine data (E37).

## Limits of this reconciliation

- **The stabilizer is not an equivalence point.** ANDES's IEEEST is structurally
  different from the internal washout-plus-lag on electrical power. The
  reconciliation is exact with the stabilizer **removed on both sides**. Results
  quoted with the stabilizer active are internal-model results and remain not
  independently reproduced.
- **The converter is not comparable at all.** ANDES's REGCA1/REECA1 are different
  devices. The reconciliation uses the static-injection configuration — which is
  where the original disagreement arose, so this is not an evasion.
- **The machine order is equivalent in band only.** GENROU is degraded with
  `xd2 = xd1`, `xq2 = xq1` and subtransient constants of 1e-4; a decade change to
  1e-3 leaves the band results bit-identical.

## A stale artifact found and corrected

`configs/ias2026/ieee39_ybus_andes.npy` differs from the internal Ybus by exactly
`−1j` at bus 4 and `−2j` at bus 5, i.e. it is a **shunt-free** network. The live
ANDES case carries `Shunt_1 b = +1.0` at bus 4 and `Shunt_2 b = +2.0` at bus 5,
and the internal Ybus includes them. The stored file is stale and does not match
the reference power flow it sits beside; the live ANDES data now written to
`results/F1/andes_live_powerflow.json` should be used instead. This did not
affect any result — the internal Ybus was always the correct one — but the stored
array should not be used as an ANDES reference.

## Gate outcome

> If they agree after harmonization: freeze the harmonized benchmark and make it
> the ONLY IEEE-39 model eligible for journal claims.

They agree. `configs/ieee39_harmonized_dynamic_model.yaml` is frozen and is now
the only eligible IEEE-39 dynamic model. TPWRS claim escalation is **not**
blocked, subject to the excitation conditionality above.

## Files

- `docs/F1_ANDES_MODEL_RECONCILIATION.md` — this document
- `results/F1_jacobian_block_comparison.csv` — the 17-item subsystem comparison
- `results/F1_eigenvalue_reconciliation.csv` — per-case band comparison
- `results/F1/F1_band_spectrum_match.csv` — band-restricted spectrum matching
- `results/F1/F1_andes_services.csv` — the service-ablation ladder
- `results/F1/F1_andes_equivalent.csv`, `F1_andes_spectra.npz` — equation-equivalent stage
- `results/F1/andes_live_powerflow.json` — live ANDES network and solution
- `configs/ieee39_harmonized_dynamic_model.yaml` — the frozen benchmark
