# E31 — independent ANDES validation: **G2 MIXED**

ANDES 2.0.0, its own interpreter, its own case data, none of this project's code.
Different machine model (GENROU, six states with damper windings, against the
four-state model here), different exciter (IEEEX1), different stabilizer
(IEEEST), different network assembly, different eigenvalue path.

Bands were declared in the experiment source **before** any lattice comparison
was inspected: base frequency within 10 % relative, base damping sign must agree,
the full portfolio must be worst in both, and the per-subset shift sign must
agree in at least 80 % of the fifteen non-empty subsets. Magnitudes were never to
be compared, because the device models differ by construction.

## What is comparable, and what is not

ANDES has no model of this project's grid-following converter. Its REGCA1 and
REECA1 are structurally different devices, so using them would compare two
converters rather than validate one. The shared configuration both tools can
express is **machine removed, power still delivered as a constant-power
injection** — exactly the negative control this project calls `static_power`, and
the configuration in which `CLAIMS.md` N4c asserts the failure persists with no
converter present. That is what was compared. The converter repair cases (RC, RB)
are reported **NOT COMPARABLE**, not compared badly.

## Result

| | this work | ANDES | agreement |
|---|---|---|---|
| base inter-area frequency | 0.6385 Hz | 0.6645 Hz | **4.07 %**, inside the band |
| base inter-area damping | −0.1265 | −0.2064 | **sign agrees**; ANDES more damped |
| median frequency error over the lattice | | | **6.2 %** |
| flagship, static injection | **+0.2624, unstable** | **−0.3323, stable** | **disagree** |
| full portfolio is the worst subset | yes | **no** | **disagree** |
| Spearman of α against subset size | −0.34 | **−0.93** | **disagree in trend** |
| per-subset shift sign agreement | | | **73.3 %**, below the 80 % band |

**In ANDES, removing these four machines makes the inter-area mode monotonically
more damped.** In this work it makes it unstable. The two implementations
disagree on the *sign* of the effect, not merely its size.

## What this does and does not overturn

It does **not** test the grid-following flagship: ANDES cannot represent it.
Nothing here bears directly on the GFL results of v2A/v2B/v2C.

It **does** contradict the static-injection configuration, which is the evidence
base for `CLAIMS.md` N4c, "the failure persists with no converter present". An
independent implementation of that same configuration finds no failure at all.

## Candidate causes, none yet eliminated

1. **Machine order.** GENROU carries damper windings; this project's machine is
   fourth order. Suppressing the damper windings in ANDES (`xd2 = xd1`,
   `xq2 = xq1`, subtransient constants to 1 ms) moves the ANDES base from
   −0.2064 to −0.1634 at 0.6475 Hz, i.e. **within 1.4 % of this work in frequency
   and much closer in damping** — but the flagship still becomes *more* damped,
   −0.1943. So the disagreement survives this correction.
2. **Excitation.** ANDES's IEEEX1 data carries `KE = −0.05`, which gives the stock
   case **six unstable real modes with a spectral abscissa of +1.03** before any
   replacement. This project rejected IEEEX1 for exactly that reason and uses a
   first-order AVR with the case `KA`/`TE`. The independent implementation is
   therefore running an excitation model this project already judged defective.
   Those six modes lie outside the inter-area band and do not enter the tracked
   comparison, but they show the ANDES dynamic dataset is not clean.
3. **Stabilizer.** IEEEST against a washout-plus-lag on electrical power.
4. **An error in this project's static-injection path**, not excluded.
5. **An error in the ANDES case rebuild**, not excluded.

## Decisive experiment, not run tonight

Write a custom ANDES model card implementing this project's device equations —
fourth-order machine, first-order AVR, power-input PSS, and the ten-state
grid-following converter — and re-run the lattice. That converts "two different
models disagree" into "two implementations of the same model agree or do not".
This is recorded in `TPWRS_REMAINING_GAPS.md` as the highest-priority remaining
gap for a Transactions submission.

## Claim actions taken

- `N4c` downgraded: its evidence is **not independently reproduced**. The
  statement "the failure persists with no converter present" now carries a
  MIXED independence status.
- `C1`/`C2` (proper subsets stable, portfolio unstable) keep their internal
  evidence but are marked **not independently confirmed** for the poster and for
  TPWRS.
- The base-case inter-area mode **is** independently confirmed, in frequency and
  in damping sign, and that may be stated.

## Files

`E31_ANDES_validation.csv`, `E31_andes_raw.csv` (32 ANDES cases, both governor
settings), `E31_ANDES_vs_internal.png`, `E31_figure_source.csv`, `cases/`
(the 16 rebuilt ANDES workbooks), `manifest.json`.
