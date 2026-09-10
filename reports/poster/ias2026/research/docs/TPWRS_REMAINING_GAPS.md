# Remaining gaps for a Transactions on Power Systems submission

Ordered by how likely a reviewer is to stop on them. Each names the decisive
experiment, not a wish.

## Gap 1 — independent implementation, unresolved (blocking)

**Status: G2 MIXED.** ANDES 2.0.0 reproduces the base inter-area mode — 0.6645 Hz
against 0.6385 Hz, 4.1 %, damping sign agreeing — and does **not** reproduce the
machine-removal ordering. Removing the four machines makes the mode monotonically
*better* damped in ANDES (−0.206 to −0.332) and unstable here (−0.126 to +0.262).
The two implementations disagree on the **sign** of the effect.

Not yet eliminated: the six-state GENROU with damper windings against the
four-state machine used here (suppressing the damper windings moves ANDES to
−0.163 at 0.6475 Hz, within 1.4 % in frequency, and the flagship still becomes
better damped); ANDES's IEEEX1 data carrying `KE = −0.05`, which gives its own
stock case six unstable real modes before any replacement; the IEEEST stabilizer
against the washout-plus-lag used here; an error on either side.

**Decisive experiment.** Write a custom ANDES model card implementing this
project's device equations — fourth-order machine, first-order AVR, power-input
PSS, and the ten-state grid-following converter — and re-run the 16-subset
lattice. That converts "two different models disagree" into "two implementations
of the same model agree or do not". Until then no claim beyond the base mode may
be called independently validated.

## Gap 2 — the machine-parameter sensitivity is unmapped

E37 shows the full pattern holds in **418 of 1000** draws over ±20 % inertia,
±10 % transient reactance, ±15 % AVR and ±20 % PSS. A reviewer will ask which
region of that space supports the result.

**Decisive experiment.** A sensitivity decomposition — Sobol indices or a fitted
response surface on the order-4 indicator over the six machine coordinates — to
identify which parameters carry the structure and over what range. The existing
1000 draws are enough to fit it; nothing new needs simulating.

## Gap 3 — the repair is characterised only at one operating point

E34 finds a converter-only solution at every declared margin **at the nominal
point**. E35 finds the *frozen* retune works at 117 of 240 held-out unstable
points. Nobody has asked whether an *adaptive* retune — re-optimised per
operating point — succeeds everywhere.

**Decisive experiment.** Run the E34 M1 optimisation at each of the 240 unstable
E35 samples and report the success rate and the required change norm as a
function of severity. Roughly 240 SLSQP runs; hours, not days.

## Gap 4 — dispatch dependence needs a proper map

E30 tested four reactive policies at one operating point and found the
interaction order to be 4, 3 and absent. That is a striking result on a single
point.

**Decisive experiment.** Repeat the four-policy 16-subset lattice across the E33
continuation grid, so the dispatch dependence is reported as a region, not an
anecdote. Cost is four times E33, which ran in 20 s on a quiet machine.

## Gap 5 — the bus-30 audit is underpowered

E38's unadjusted evidence is strong (12 of 12 genuine inter-area cores, matched
pairs 41–0, McNemar p = 9.1e-13) and the adjusted odds ratio is 4.09 with an
interval of [0.14, 124] on ten unstable cases.

**Decisive experiment.** Repeat the audit over the per-sample lattices already
stored in E35, which contain 240 unstable operating points rather than 10. The
data exist; only the analysis is missing.

## Gap 6 — extended precision was never available

`numpy.longdouble` has a 53-bit mantissa on this platform and no package could be
installed. E40 substitutes three solve algorithms, three perturbation amplitudes
and three grid densities, and shows floating point contributes at most 2.8e-10
against a smallest margin of 0.00488 — a headroom factor of 1.7e7 — and that the
**frequency grid dominates by eight orders of magnitude**.

**Decisive experiment.** An `mpmath` recomputation of `m₄` at five near-boundary
points, once a package may be installed. Low risk; the perturbation ensembles
already bound the answer.

## Gap 7 — the closure metric is grid-limited near the boundary

The frozen 61-point metric is a conservative upper bound whose looseness
concentrates exactly where the margin is small: refined, the smallest audited
margin falls from 0.00488 to 0.00116. The bias runs *against* every claim made
from it, so nothing is overstated, but a paper should present the refined metric
as primary with the frozen one as the preregistered check.

**Decisive experiment.** Preregister the locally refined band minimum as the
primary metric in a successor protocol and re-run E33 and E35 on fresh seeds.

## Gap 8 — one benchmark, one topology

Everything here is IEEE-39. Nothing licenses a claim about other systems.

**Decisive experiment.** Repeat the flagship search on a second benchmark — IEEE
118 or a synthetic Texas case — and report whether a minimal incompatible core of
order 4 exists at all.

## Not gaps

- **Numerical trustworthiness of the closure.** E40 settles it within the
  precision available.
- **Port reduction fidelity.** E41 settles it at 1e-10 across the lattice.
- **Linear-to-nonlinear consistency.** E32 settles it at 0.003 Hz.
- **Controller-cause.** E36's 1500 of 1500 settles it as well as a Monte Carlo can.
