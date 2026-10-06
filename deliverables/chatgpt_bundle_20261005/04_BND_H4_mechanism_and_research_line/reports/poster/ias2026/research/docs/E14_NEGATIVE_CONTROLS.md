# E14 — Negative controls

394 cases, 43 s. Artifacts: `results/tables/E14_negative_controls_*.csv`,
manifest `results/manifests/E14_negative_controls.json`.

This experiment can only falsify. Its result changes the flagship selection.

## Hard gate — PASSED

The method is not an "everything is dangerous" detector.

- Order-4 flag rate: **10 of 126 = 7.9 %**.
- **N4**: 25 stable four-replacement portfolios matched to the failing ones —
  replaced MW 4665–5025 against a failing range of 4144–4985, identical minimum
  SCR (0.907), gSCR 0.833–0.897 — are **all stable**, worst spectral abscissa
  `-0.0677`. Matched size, matched MW, matched conventional strength, opposite
  outcome.

## N5 — the controller-mediation test, and it splits the two families

Each machine is replaced by an ideal injection at **the same P and Q**, with no
dynamics. Two bracketing static references are used: constant power (harshest —
its negative incremental impedance is destabilizing on its own) and constant
impedance (most benign). A constant-current variant was implemented and then
removed: with the injection angle following the bus voltage it leaves the
voltage magnitude undetermined and the equilibrium degenerate.

Median spectral abscissa:

| failure family | GFL | static P/Q | static impedance | controller-mediated |
|---|---|---|---|---|
| genuine inter-area (0.42–0.57 Hz), 9 cores | +0.1447 | **+0.2624** | +0.1433 | **0 of 9** |
| mode-switching / PLL (f ≈ 0), 30 cores | **+322.19** | **−0.1073** | +0.2913 | **23 of 30** |

**The inter-area order-4 family is NOT controller-mediated.** Every one of the
nine cores is unstable when the four machines are replaced by injections with no
dynamics at all, under both references, and the converters are consistently
*less* destabilizing than the static reference (by 0.016 to 0.188, mean 0.091).
The mechanism is the removal of four synchronous machines — their inertia,
damping torque, AVR and PSS — from the inter-area mode.

**The PLL family IS controller-mediated, decisively.** With the same P and Q and
no dynamics the system is stable at `-0.197` to `-0.237`, which is a *larger*
margin than the base case at `-0.126`. Adding the converters takes it to `+185`
to `+1276`. These 23 cores are of size 5 (13) and size 6 (10). The size-4 PLL
case `33+35+36+38` is not among them: its static replacement is also unstable at
`+0.300`, so the smallest PLL core is contaminated by the inter-area effect.

## N1 — the inter-area mechanism is discontinuous in the replacement fraction

Spectral abscissa of the nine genuine cores against `rho`:

| rho | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|
| range over the nine cores | −0.144 to −0.166 | −0.175 to −0.229 | −0.190 to −0.237 | **+0.027 to +0.389** |

The mechanism does not emerge continuously. Damping **improves monotonically**
up to `rho = 0.75` and then jumps unstable at `rho = 1.00`, where the machine is
removed entirely rather than merely shrunk. It is a binary function of whether a
synchronous machine remains at the bus, not a function of converter penetration.

This is an independent confirmation of N5 and it rules out a
penetration-driven-interaction reading of this family.

## N2 — synchronous damping does not repair it

The source case has `D = 0` on all ten machines, so a multiplicative "increase
damping" control is vacuous. Damping was added absolutely, `D ∈ {0, 1, 2, 5}` pu
on machine base, inertia untouched.

Effect on the spectral abscissa at `D = 5`: between `-0.026` and `+0.029`. For
four of the nine cores it makes matters **worse**. Damping torque is not the
mechanism.

## N3 — controller ablation: not controller-caused, but controller-repairable

Median spectral abscissa by ablation:

| ablation | inter-area family | PLL family |
|---|---|---|
| outer P loop ÷4 | **−0.185** | +389.6 |
| outer Q loop ×4 | +0.116 | +519.0 |
| measurement fast | +0.141 | +519.5 |
| nominal | +0.145 | +413.2 |
| PLL slow (wn 12) | +0.175 | **+190.3** |
| outer P loop ×4 | +0.284 | +498.2 |
| outer Q loop ÷4 | +0.319 | +382.8 |

Per core, weakening the active-power outer loop by a factor of four **stabilizes
9 of 9** genuine cores, with margins from `-0.032` to `-0.235`, and the critical
mode stays in the synchronous family in all nine.

So the two statements are both true and must not be conflated:

- **not controller-caused** — the failure exists without any converter (N5, N1);
- **controller-repairable** — a converter outer-loop retune removes it while all
  four PV replacements stay connected (N3).

For the PLL family, slowing the PLL from `wn = 37.4` to `12 rad/s` more than
halves the divergence but never stabilizes it, consistent with the earlier
finding that it is not a bandwidth problem.

## Consequences for the flagship

The chain this phase set out to test —

    all subsets safe -> irreducible order-4 -> controller-mediated port core
    -> provenance -> targeted repair

**cannot be run on the inter-area family as specified.** N5 refutes the
controller-mediation link for all nine cores, and N1 refutes the continuity that
a penetration story would require. Continuing to E15/E16 on that family under a
controller-mediated framing would be building on a link the controls just broke.

Two viable tracks remain, and they answer different questions:

**Track A — inter-area, reframed.** Keep `30+33+35+37` as flagship, drop every
mediation claim, and state it as: four full machine retirements that are
pairwise and triple-wise safe destabilize the inter-area mode through a
genuinely irreducible fourth-order term, and a converter outer-loop retune
repairs it without undoing any replacement. Supports C1, C2, C6. Does not
support a controller-mediated core claim. The Moebius order decomposition is
valid here because the mode family is preserved across the lattice (min MAC
0.967).

**Track B — PLL, where mediation is established.** Take a flagship from the 23
controller-mediated cores, best static margins being `31+32+35+36+38` (static
`-0.237`, GFL `+908`) and `32+35+36+37+38` (static `-0.232`, GFL `+1065`).
Here N5 gives clean evidence of controller mediation. The cost is that the
Moebius order decomposition is invalid for these cores — the critical mode
changes family across the lattice, which is why E12 classified them
`C_MODE_SWITCHING` — so the order claim must be carried by **winding
provenance**, which does not require mode identity, rather than by Moebius.

Running both is defensible and cheap; the port-space machinery is shared. Track A
gives the order-4 Moebius waterfall figure, Track B gives the controller-mediated
provenance and repair figure. Neither alone supports the full original chain.

## N6 — the matched stable controls do not hide the same branch

The flagship inter-area branch of `30+33+35+37` was tracked into the matched
stable portfolios by MAC, computed on the machines surviving in **both**
portfolios (4 or 5 machines; the naive comparison fails because the two
portfolios retire different machines and therefore have different state vectors).

15 of the 25 matched controls share at least four machines with the reference.
MAC to the reference branch: minimum 0.941, median 0.992 — the same physical
mode is being followed. Its real part in those controls is at most `-0.0677`,
with damping at least `2.17 %`.

So the matched controls are not near-misses in disguise: the branch that goes
unstable in the failing quadruples stays clearly damped in the matched stable
ones. The discrimination in N4 is about this mode, not about a different one.
