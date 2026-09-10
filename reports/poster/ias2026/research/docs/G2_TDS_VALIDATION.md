# Gate 2 — Nonlinear time-domain validation (phasor domain, NOT EMT)

**This is not an electromagnetic-transient study.** It is a nonlinear
phasor-domain DAE simulation. The network is its algebraic current balance and
the devices are their differential equations, integrated without
linearization:

- BDF integration;
- Newton solve of the network at every right-hand-side call;
- central-difference Jacobian.

No switching, no line dynamics, no converter PWM, no controller limits. The
question is whether the nonlinear model does what the eigenvalues say:
frequency, growth or decay, stability outcome, and which side of a threshold a
case is on.

Code: `experiments/G2_tds.py`. Data: `results/G2/` (one trace per run). The
first full run, with an integrator defect, is kept in `results/G2/first_run/`
(see §4).

## 1. Design

Two disturbances. Both leave the equilibrium unchanged (no governors):

- **D1**: a `+1e-4` pu rotor-speed kick on a machine present in every case
  (IEEE-39 bus 31, Kundur machine 1);
- **D2**: a `+2 %` active-load pulse lasting 0.2 s (IEEE-39 bus 20, Kundur
  bus 7).

Integration settings:

- Horizon `clip(3/|alpha|, 30, 120)` s.
- The observable is a machine-speed difference. For a real mode it is the
  state with the largest participation.
- Frequency and growth come from a matrix-pencil ringdown fit in the
  small-signal window (spread up to 1e-3 pu), plus an FFT and an envelope
  slope.
- The outcome is model-free: DIVERGED, UNSTABLE (the envelope of the second
  half grows) or STABLE.

16 cases x 2 disturbances = 32 runs.

## 2. Results

`alpha` is the growth rate in s^-1 and `f` the frequency in Hz. The TDS column
gives the pencil exponent nearest the predicted mode, as the D1 / D2 range.

| case | prediction (`alpha`, `f`) | TDS (`alpha`, `f`, D1 / D2) | outcome D1 / D2 |
|---|---|---|---|
| **IEEE-39** | | | |
| SAFE, P_inf flagship (`H = EMPTY`) | −0.1554, 0.7276 | −0.1555 / −0.1554, 0.7276 | STABLE / STABLE |
| `kappa = 4`, P4 flagship `30+33+35+37` | **+0.1270**, 0.6223 | +0.1269 / +0.1274, 0.6223 | UNSTABLE / UNSTABLE |
| P4 proper subset `30+33+35` | −0.2047, 0.9160 | −0.2053 / −0.2050, 0.9161 | STABLE / STABLE |
| `kappa = 2`, P2 pair `30+33` | **+0.0717**, 0.6421 | +0.0717 / +0.0717, 0.6421 | UNSTABLE / UNSTABLE |
| P2 single `30` | −0.0909, 0.6248 | −0.0912 / −0.0913, 0.6248 | STABLE / STABLE |
| P2 single `33` | −0.0241, 0.6402 | −0.0242 / −0.0241, 0.6402 | STABLE / STABLE |
| tongue line `k = 1.30`, `g = 0.020` (before the gap) | **+0.0784**, 0.5727 | +0.0790 / +0.0786, 0.5726 | UNSTABLE / UNSTABLE |
| `g = 0.042` (gap) | −0.0082, 0.6242 | −0.0082 / −0.0082, 0.6242 | STABLE / STABLE |
| `g = 0.102` (tongue interior) | **+0.0183**, 0.6808 | +0.0183 / +0.0183, 0.6808 | UNSTABLE / UNSTABLE |
| `g = 0.173` (after the tongue) | −0.0082, 0.7013 | −0.0083 / −0.0083, 0.7014 | STABLE / STABLE |
| condenser `D = 2` at 2.0 % (RHP threshold 2.48 %, G1) | **+0.0264**, 0.6209 | +0.0264 / +0.0268, 0.6209 | UNSTABLE / UNSTABLE |
| condenser `D = 2` at 3.0 % | inter-area −0.0295, 0.6213 (rightmost: real −0.0017) | −0.0294, 0.6213 (real: −0.0011 / −0.0001) | STABLE / STABLE |
| undamped classical condenser, P_inf, `37` (G1 swing mode) | **+0.0345, 8.6502** | +0.0341 / +0.0349, 8.6503 | UNSTABLE / UNSTABLE |
| **Kundur** | | | |
| oscillatory boundary, replacement `2`, `k = 1.25`, `g = 0.08` | −0.0850, 0.6072 | −0.0852 / −0.0847, 0.6072 | STABLE / STABLE |
| same, `g = 0.11` (boundary at `g = 0.0949`) | **+0.0548**, 0.6289 | +0.0548 / +0.0548, 0.6289 | DIVERGED* / UNSTABLE |
| aperiodic, replacements `2+3`, `k = 0.75`, `g = 0` | **+0.685**, real | D1: +0.77 (real, log-slope of `efd_sg1`, 1.9–5.7 s) | DIVERGED / DIVERGED** |

\* The oscillation grows at the predicted rate for 48 s. The integration then
ends (step-size failure at 50.6 s) once the speed spread reaches `4e-3` pu.
This is recorded as DIVERGED and is not interpreted further.

\** D1: exponential departure along the predicted real mode, which is an
excitation/flux mode: `efd`, `eq'` and the AVR lead-lag of machines 1 and 4
carry 73 % of the participation.
The network solution is lost at 6.4 s. D2: the 2 % load pulse is a large
disturbance at this operating point. The network solution is lost within
0.1 s, during the pulse, so D2 gives no modal comparison for this case.

**Summary**

| check | result |
|---|---|
| stability verdict (TDS outcome vs `N_RHP > 0`) | **32 / 32** |
| oscillatory cases, pencil vs eigenvalue (28 runs) | max `|df| = 2.5e-4` Hz (median `1.0e-5`); max `|d alpha| = 6.0e-4` s^-1 (median `9.6e-5`), at most 1.8 % relative; growth sign 28/28 |
| aperiodic case (D1) | sign and order of magnitude agree; rate +0.77 against +0.685 (+12 %, fitted up to where the deviation reaches 0.05, which includes nonlinear acceleration) |
| threshold behaviour, condenser rating | `alpha` of the inter-area mode interpolated between 2 % (+0.0264) and 3 % (−0.0294) is zero at **2.47 %**; bisected RHP threshold (G1) **2.48 %** |
| threshold behaviour, tongue | U → S → U → S along `k = 1.30` at `g = 0.020, 0.042, 0.102, 0.173`, as predicted |
| threshold behaviour, Kundur | S at `g = 0.08`, U at `g = 0.11`; linear interpolation gives 0.098 against the located boundary 0.0949 (`alpha(g)` is not linear) |
| `kappa` structure | `kappa = 4`: flagship grows, triple decays. `kappa = 2`: pair grows, both singles decay. |
| G1 condenser finding | the undamped condenser's own 8.65 Hz swing mode grows in the nonlinear model at the predicted rate |

One pencil artefact: in K4-proper D1 the component dominating the end of the
window is a slow 0-Hz drift (`+0.027`). The exponent at the predicted mode is
−0.2053. The outcome rule (envelope) gives STABLE, as predicted.

## 3. Verdict

The nonlinear phasor-domain model reproduces, for every declared case:

- the frequency and growth rate of the critical mode, to within `2.5e-4` Hz
  and `6e-4` s^-1 for oscillatory modes;
- the stability outcome, 32 of 32;
- which side of the condenser, tongue and Kundur policy thresholds each case
  is on.

The small-signal `H_RHP` labels of the validated cases are therefore
confirmed in the time domain. Nothing here is a transient-stability or EMT
result. Large disturbances, limits and fast electromagnetic dynamics are
outside the model.

## 4. Integrator record (FAILED_EXPERIMENTS F18)

Before the full run:

- The sustained load step (E25/E32) was replaced by the kick and the pulse.
  Without governors a sustained step has no equilibrium, and the frequency
  drifts out of the small-signal regime.
- The fixed equilibrium Jacobian was replaced, because BDF's Newton iteration
  stalled on PLL transients.

The first full run used scipy's adaptive finite-difference Jacobian:

- On the identically zero columns of the frozen condenser EMFs, it grew its
  step until the network equations had no solution. The NaN Jacobian was
  reported as ALGEBRAIC_FAILURE in the two condenser cases and the swing case.
- The final run uses a fixed-step central difference, and all 32 runs were
  repeated.
- The first run's valid cases agree with the final run to about `1e-3`.
