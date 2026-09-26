# Gate 5 — Novelty statement, tested against the evidence

## Not claimed

None of the following is new, and none is claimed:

| item | why it is not new |
|---|---|
| generalized Nyquist | classical multivariable stability test. F10: applied to every principal sub-loop it recovers `H` exactly because it *is* the per-subset computation (O78, N23). |
| the `-1` crossing, return difference, impedance/port criterion | classical; `det(I + M_S)` is the return difference of the portfolio loop (N23) |
| hypergraphs, clutters, minimal-element families | standard combinatorics. `H` is the set of minimal elements of the unsafe family (F2C). |
| non-monotone damping, or voltage support that destabilizes | known in the converter-interaction and PSS literature. We report it (N22) and do not claim it. |
| eigenvalue maps over parameter planes, participation factors, mode tracking | standard small-signal practice |
| the IEEE-39, Kundur and 68-bus models | documented benchmarks, used as published (68-bus reproduces Table 4 to `5e-4` Hz and `5e-4` percentage points) |

## Candidate contribution, and the test of each part

> *A policy-dependent minimal spectral incompatibility structure over a
> combinatorial family of SG-to-IBR replacement portfolios, with
> piecewise-constant regions, witness-coalition changes, reduced
> port-boundary localization, and physically actionable transitions.*

Each qualifier was tested on the primary object `H_RHP` (F2E). The band object
counts only where it coincides with `H_RHP`.

| qualifier | test | evidence | verdict |
|---|---|---|---|
| **minimal spectral incompatibility structure** | an antichain of minimal unsafe portfolios whose absence certifies safety | F2C Props 1–2; F2E Prop 1; "no hyperedge contained implies stable" holds by construction; heredity fails, so the object is not a simplicial complex | **supported** (definition plus proofs) |
| **over a combinatorial family of SG-to-IBR replacement portfolios** | the object is evaluated on the power set, not only on the flagship | 16 (IEEE-39), 8 (Kundur), 16 plus 4 096 (68-bus) portfolios; F10: flagship-only Nyquist misses the coalition structure | **supported** |
| **policy-dependent** | `H_RHP` changes with the converter reactive policy at fixed network and dispatch | IEEE-39: F7A 30, F7B 36, F7C 16 distinct `H_RHP` with `H_RHP = H_IA` at 100 / 96.5 / 100 %. Kundur: `H_RHP` changes on 61/61 and 36/36 pure-policy lines but is never empty. **68-bus: `H_RHP = EMPTY` at all 961 nodes of the preregistered 4-candidate map (NOT REPRODUCED).** The preregistered 12-plant check does show policy dependence: `kappa_RHP = 6, 9, 11` and 259, 1, 2 hyperedges at `g = 0, 0.25, 1`, with a change of mechanism from inter-area to converter PLL. That is three points, not a map. | **supported on IEEE-39 and Kundur; on the 68-bus system only at high penetration and only at three points; system-dependent** |
| **piecewise-constant regions** | local constancy away from `∂Gamma_RHP`; regions bounded only by spectral crossings | F2C Prop 3 / F2E §3; 795 RHP boundaries located on IEEE-39 and Kundur; 0/1 153 spot-check errors (F7); 60/60 (68-bus) | **supported** |
| **witness-coalition changes** | the minimal coalition itself changes, not only `kappa` | F7 contraction `{30,33,35,37} -> {30,33,35} -> {30,33}`, identical under RHP on F7A/F7C; up to 15 hypergraphs share one `kappa`; Kundur `2\|3+4` vs `2+3\|2+4\|3+4` at `kappa` 1 vs 2 | **supported** (IEEE-39, Kundur) |
| **reduced port-boundary localization** | every oscillatory boundary is a zero of the reduced return difference | 29 851/29 851 (F7 band), 1 511/1 511 (F8C), 601/601 (G1, RHP). **Aperiodic crossings through the origin are not port-certifiable** (reference zero at `s = 0`), and there are 194 of them on Kundur. The per-subset test is classical GN (N23). The reduction buys speed only at a common operating point (N24). | **supported for oscillatory boundaries only**. The novelty is the use of the classical test as a boundary locator for the hypergraph, not the test itself. |
| **physically actionable transitions** | transitions moved by physical interventions, in both directions | F8 under RHP (G1): a condenser with damped swing mode (`D = 2` or flux dynamics) empties `H_RHP`; thresholds 0.24 % (P_fold) to 45 % (P1) of the retired rating; AVR speed orders the regions; surviving-fleet inertia destabilizes; the undamped condenser is **not** safe. G2: the nonlinear phasor-domain model confirms the frequency, growth sign and threshold side of the key cases | **supported on IEEE-39**. Not tested on Kundur or the 68-bus system. |

## Final wording

> We define and compute, for a combinatorial family of SG-to-IBR replacement
> portfolios, the minimal spectral incompatibility hypergraph `H_RHP`: the
> minimal portfolios that place an eigenvalue in the open right half-plane.
> On IEEE-39 and on the Kundur two-area system it depends on the converters'
> reactive policy.
>
> - It is piecewise constant, bounded only by spectral crossings, and it
>   reorganizes through changes of the witness coalition that the scalar
>   order `kappa` does not show.
> - Every oscillatory boundary is localized by the classical return difference
>   of a reduced port model.
> - On IEEE-39, physical services (electromagnetic presence with a damped
>   swing mode, the excitation speed of the surviving fleet) move the
>   transitions in both directions.
>
> The per-subset test is generalized Nyquist, and it is not new.
>
> On the documented 68-bus NETS–NYPS system, every portfolio of the four
> preregistered candidates is stable throughout the tested policy and
> excitation window. Incompatibilities appear only when at least six of the
> twelve physical plants are replaced, and there too they depend on the
> policy. The structure is therefore system-dependent, not universal.

Mechanism-resolved projections (`H_IA`) are reported as explanations of the
declared inter-area mechanism, never as safety statements.
