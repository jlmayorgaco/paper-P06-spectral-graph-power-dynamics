# Failed and downgraded in the final validation run

Kept because they constrain what may be claimed. Run
`outputs/ias2026/final_validation_overnight_20260910T003225/`.

## F10 — the irreducible fourth-order character does not survive every dispatch (E30)

Expected the order-4 structure to be a property of the replacement. It is a
property of the replacement **and its reactive dispatch**.

| policy | class | effective order | flagship |
|---|---|---|---|
| Q-matched (mechanism isolation) | A | **4** | +0.1447, unstable |
| unity power factor | B | **3** | +0.4436, unstable and worse; the triple 30+33+35 is already unstable at +0.1407 |
| voltage regulation, 0.95 pf capability | D | — | **−0.0953, stable**, zero band RHP modes |
| voltage regulation, apparent-power capability | D | — | **−0.0953, stable** |

Under unity power factor "every proper subset is safe" is **false**. Under
voltage regulation there is no instability to explain. Capability constants were
frozen before any policy case was solved.

**Consequence.** The fourth-order claim may never be stated without its dispatch.
Recorded as O48.

## F11 — the independent implementation does not reproduce the ordering (E31)

ANDES 2.0.0 reproduces the base inter-area mode to 4.1 % in frequency with the
damping sign agreeing, and reproduces **none** of the machine-removal ordering:
per-subset shift-sign agreement 73.3 % against a declared 80 % band, the full
portfolio is not the worst subset in ANDES, and the flagship is **stable at
−0.3323** there against **+0.2624** here.

This bears on the static-injection configuration, which is the evidence base for
`CLAIMS.md` N4c, "the failure persists with no converter present". It does not
test the grid-following flagship, which ANDES cannot represent.

**Consequence.** N4c downgraded; C1 and C2 marked not independently confirmed.
Recorded as N13.

## F12 — the fixed converter retune is not a robust repair (E35)

Expected the E21 retune to generalise. Across 240 held-out unstable operating
points it restores stability in **117**, 0.4875 [0.4227, 0.5526], and its
authority is bounded in an interpretable way:

| severity | RC succeeds |
|---|---|
| `α_IA` ≤ 0.122 | **80 of 80** |
| 0.122 < `α_IA` ≤ 0.242 | 37 of 80 |
| `α_IA` > 0.242 | **0 of 80** |

The 25 % synchronous condenser succeeds in **240 of 240**.

**Consequence.** "A converter retune repairs stability while retaining all PV" is
withdrawn as an unconditional claim. Recorded as N14.

## F13 — the order-4 structure is not robust to machine data (E37)

Over ±20 % inertia, ±10 % transient reactance, ±15 % AVR and ±20 % PSS, the full
pattern holds in **418 of 1000** draws, and "all proper subsets stable" in 721.
Not explained by a degenerate base: the base case is unstable in only 2 of 1000
draws and conditioning changes nothing.

The contrast with E36 — 1500 of 1500 under controller uncertainty — is itself the
result: the effect is insensitive to the converters and sensitive to the machines,
which supports the physical reading and limits the generality.

**Consequence.** The order-4 result is conditional on this machine data.

## F14 — the bus-30 collective-enabler hypothesis is not supported (E38)

Unadjusted, the evidence looked strong: bus 30 is in **12 of 12** genuine
inter-area minimal cores and is the only bus at 100 %; Fisher odds ratio 13.2,
p = 0.0050; matched pairs **41 wins, 0 losses**, McNemar exact p = 9.1e-13;
permutation p = 0.0026. Adjusted for megawatts, inertia and short-circuit
strength the odds ratio is **4.09 [0.14, 124.2], p = 0.418**.

The preregistered rule — do not claim it if bus 30 loses significance after
controlling for megawatts and inertia — is not met. That the adjustment is
underpowered (ten unstable cases, five predictors, the restricted fit failing to
converge) does not rescue it.

**Consequence.** Descriptive statement only. Recorded as N15.

## F15 — "conventional baselines are anti-predictive" was an overstatement (E39)

Removed inertia ranks instability at AUC 0.78–0.87 and replaced megawatts at
0.77–0.86. Only the short-circuit family is weak, 0.59–0.72 in the conventional
direction. The earlier framing was too broad.

What survives, and more strongly: under the declared threshold the lower-order
reconstructions produce **0 true positives and 0 false positives at every size**.
They do not rank badly — they declare every portfolio stable, including the 139
that are not.

**Consequence.** N16 corrects the baseline claim; N17 records the stronger
threshold result.

## Three code defects found and fixed during the run

**E30 classifier.** Decided "mechanism disappears" on whether damping worsened
relative to base, which labelled the voltage-regulating policies B although their
flagship is stable with no right-half-plane mode in the band. Now decided on the
absence of an instability. No threshold or model changed.

**E34 target test.** Judged an optimizer that converges *onto* its constraint at a
`1e-9` tolerance, so a successful repair read as a failure. At `1e-6` the
converter-only strategy meets all four targets. This inverted the headline of the
experiment, which is why it was worth re-running.

**E40 perturbation ensemble.** Read the closure eigenvalue off the *unperturbed*
operator, so the reported eigenvalue variation was identically zero — a metric
measuring nothing. Corrected; the margin statistics were unaffected.

## Not failures

- **H5B** remains NON-EXECUTABLE from v2B. Nothing tonight re-scores it.
- **v1** remains frozen and refuted; **v2A** and **v2B** remain frozen.
- The E34 finding that a marginal repair has a *smaller* closure margin than the
  unrepaired flagship is **not** a failure of the invariant: `m₄` measures
  distance to the boundary in either direction, exactly as E33 established.
