# Symmetry-deflated zero-frequency port closure

Status:
- The construction and the determinant identity are proved (algebra). The
  ingredients are classical: Brauer's rank-one eigenvalue theorem and the Schur
  determinant formula.
- The frozen definition is `configs/port_relocation_v1.yaml`, including
  amendment v1.1, committed before the holdout ran.
- The validation is the holdout (FC02) in
  `results/zero_frequency_port_validation.csv`.
- The BC01 result 194/194 on the G1 Kundur crossings is **method development**:
  the relocation was designed on those same crossings. It is never cited as
  validation.

## 1. Why the plain port is blind at s = 0

For the phasor DAE `x' = f(x, z)`, `0 = g(x, z)`, the bus-voltage port is the
Schur complement

    T(s) = g_z − g_x (sI − f_x)^{-1} f_z,
    det T(s) · det(sI − f_x) = det(g_z) · det(sI − A),   A = f_x − f_z g_z^{-1} g_x.

The exact center subspace `C = span{R_x, w}` (`TRANSVERSE_STABILITY_QUOTIENT.md`
§2) gives `det(sI − A) = s^2 det(sI − A_perp)`. As a result:

- `det T(0)` vanishes structurally at **every** operating point, for every
  portfolio and every policy;
- its sign carries no information about a physical real eigenvalue crossing the
  origin;
- this is the failure BC01 observed on the Kundur Q/V-integrator crossings.

## 2. Construction (frozen, v1)

Let `U = R_x / ||R_x||`, with `beta, beta2 > 0` (frozen at 1 s^-1, consistency pair
2 and 0.5).

1. `f_x# = f_x − beta U U^T`. Then `A# = A − beta U U^T`, because the update
   enters `A` unchanged.
2. `x0 = w + a R_x`, with `a = omega_B / beta − (U^T w) / ||R_x||`.
3. `f_x## = f_x# − beta2 x0 x0^T / ||x0||^2`, so
   `A## = A# − beta2 x0 x0^T / ||x0||^2`.

Port at zero:

    T##(0) = g_z − g_x (f_x##)^{-1} f_z.

**Lemma 1 (exact relocation).** `sigma(A##) = sigma(A_perp) ∪ {−beta, −beta2}`,
with multiplicities.

*Proof.*
- Brauer's theorem: if `A u = lambda u`, then `A − u v^T` has spectrum
  `{lambda − v^T u} ∪ (sigma(A) \ {lambda})`, counting multiplicity.
- Step 1: `u = R_x` (`A R_x = 0`) and `v = beta R_x / ||R_x||^2` move one copy of
  0 to `−beta`. The second copy, from the Jordan chain, stays at 0.
- Its eigenvector for `A#` is exactly `x0`:
  `A# x0 = omega_B R_x − beta U (U^T w + a ||R_x||) = 0` for the stated `a`.
- Step 2: Brauer again with `u = x0` moves it to `−beta2`.
- No other eigenvalue moves, so the remaining spectrum is
  `sigma(A) \ {0, 0} = sigma(A_perp)`. ∎

**Lemma 2 (port identity).** Whenever `g_z` and `f_x##` are nonsingular,

    det T##(0) = det(g_z) · det(A##) / det(f_x##)

(Schur determinant formula for the block Jacobian).

**Proposition 3 (what the sign detects).** Along a continuous parameter path on
which `det g_z` and `det f_x##` keep their signs, `sign det T##(0)` changes
between two points **iff** the number of positive real eigenvalues of `A_perp`
changes parity.

*Proof.*
- Complex pairs contribute `|lambda|^2 > 0` to `det A_perp`.
- The two relocated eigenvalues contribute `beta · beta2 > 0` to `det A##`.
- Hence `sign det A## = (−1)^{#negative real eig of A_perp}`.
- The number of real eigenvalues has the parity of `dim A_perp`, so this equals
  `(−1)^{dim A_perp} (−1)^{#positive real eig}`.
- Lemma 2 and the fixed signs of the device factors complete the proof. ∎

The frozen detection rule is this proposition applied to the preregistered
brackets. It adds two requirements: no device flip (`sign det(−f_x##)` equal at
both ends) and the same verdict for both beta pairs.

**Scope, stated as limitations.**
1. It detects crossings **through s = 0** (aperiodic, real). It is blind by
   construction to oscillatory (Hopf) crossings. Those are the negative controls
   of the holdout, and the port must *not* flip on them.
2. It detects **parity**. Two simultaneous real crossings, or two positive real
   eigenvalues coalescing into a complex pair (no eigenvalue passes through 0),
   produce no flip. For a coalescence this is correct: nothing crossed the origin.
3. It needs the verified partner `w`, i.e. `D = 0` and no governor. For
   `C = span{R_x}` (the damped condensers, the governed model) only step 1 is
   needed. That variant is **not validated** here.
4. It says nothing about the frequency-domain part of the whole-RHP count. It
   closes only the real-axis point where the structural symmetry blinds the port.

## 3. Holdout validation (FC02)

The Kundur slices were never used to design the relocation:
- KH1: `t = 1.6`, five `k` values;
- KH2: `t = 1`, six `k` midpoints;
- KH3: `t = 0.7`, three `k` midpoints.

Each slice covers all 8 subsets of {2, 3, 4}, with a real scan `g in [0, 0.003]`
(61 points). The ground truth is the eigenvalue count of `A_perp`, which is
independent of the port.

| source | real crossings | detected | missed | negative controls | false positives | beta-inconsistent | max identity residual |
|---|---|---|---|---|---|---|---|
| Kundur holdout | 29 | 28 | 1 | 445 (109 oscillatory, 336 no-crossing) | 0 | 0 | 8.8e-12 |
| synthetic (200 draws, seed 20260921) | 200 | 200 | 0 | 200 | 0 | 0 | 6.0e-8 |
| IEEE-68 | 0 eligible | – | – | – | – | – | – |

**The single miss** is at KH1, `t = 1.6`, `k = 1.45`, subset 2+4, bracket
`g in [1.5e-4, 2.0e-4]`:
- the real-RHP count goes 3 → 1 and the complex-RHP count goes 0 → 2;
- two positive real eigenvalues coalesce into a complex pair, and no eigenvalue
  passes through 0;
- by Proposition 3, a flip there would have been an error.

Under the preregistered rule it still counts as a miss (28/29). The explanation
is reported beside it, not instead of it. Restricted to crossings through the
origin, the holdout result is 28/28.

Conditioning:
- the smallest relocated eigenvalue at the fixed bracket endpoints is 5.8e-5;
- amendment v1.1 replaced evaluation at 40-step bisection endpoints, where the
  crossing eigenvalue is about 1e-11 and no test is decidable, including the
  ground truth.

**Benchmark relevance.**
- The aperiodic (Q/V-integrator) crossings are a Kundur phenomenon.
- The IEEE-39 incompatibility map is oscillatory (the 0.64–0.71 Hz family).
- IEEE-68 has no real transverse crossing in the frozen map.
- The port closure is therefore a **method-level** result, validated on a
  held-out Kundur set and on synthetic descriptor systems. It is not an IEEE-39
  result.

## 4. Allowed wording

> "Because the phasor DAE is rotation-invariant and, without primary frequency
> control, admits an exact common-frequency drift, the bus-voltage port has a
> structural double zero at s = 0. Relocating exactly these two center
> eigenvalues by two rank-one updates restores the zero-frequency determinant
> test: its sign changes iff the parity of the positive real transverse
> eigenvalues changes. On a preregistered Kundur holdout it detected 28 of 29
> real-count changes (the miss is a real-pair coalescence, not a zero crossing)
> with 0 false positives on 445 negative controls."

Not allowed:
- "validated on 194 crossings";
- "detects all instabilities";
- "a new stability criterion". The ingredients (Brauer, Schur, the
  generalized-Nyquist zero test) are classical. The contribution is their exact
  application to the structural center subspace of the phasor DAE.
