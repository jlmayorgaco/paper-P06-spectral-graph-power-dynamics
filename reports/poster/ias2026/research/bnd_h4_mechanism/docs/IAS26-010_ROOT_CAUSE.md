# IAS26-010 — root cause recorded before code changes

**Recorded:** 2026-09-26, before modifying implementation code.  
**Frozen case:** H4/P4, GFL11, matched-Q policy.  
**Reproduction:** `results/20260926T100829_d0fecb32_ias26_010_reproduce`  
**Source commit:** `d0fecb3264aeb855ab6700fc6ef66120d4b6c33c`; reproduction manifest reports a clean tree.
**Pre-correction runner SHA-256:** `517a278fb81b8118ff21be9a4d641e9a7048b0a304b0ecc4b73fd11f4e46da63`.

## Reproduced first failure

The frozen runner reproduces `base+deltaY` as the first degrading stage:

- absolute residual: `7.188392772871091e-05`;
- relative residual: `2.7385662723109105e-08`;
- denominator `||T_H-T_0||_F`: `36.756566014466635`;
- frozen critical root: `0.12700646777028055 + 3.909898476444902i`.

The Ared eigenpair, descriptor pencil, raw network Schur, and candidate-port
Schur stages pass before this mismatch.  The candidate-port Schur relative
residual is `1.3083705992724766e-11` against its declared `1e-8` stage gate.

## First divergent block and mechanism

`build_action_space` independently calls `build_port_operator` at
`base_case.equilibrium.(x,z)` and `flagship_case.equilibrium.(x,z)`.  The
matched-Q contract makes their power-flow solution identical, but the two
different-sized DAE equilibrium solves return slightly different numerical
endpoints.  Consequently, the operators are not evaluated at one common
linearization point even though the rank-eight replacement identity assumes
that every non-replaced block is common.

The exact power-flow voltage arrays, network Ybus, bus ordering/index, load
model (`power`), and dispatch are identical.  The solved DAE voltage endpoints
are not: maximum distance from the common PF voltage is `5.626146037898512e-08`
for the base and `1.7031048696991158e-07` for H4; their maximum mutual voltage
difference is `2.2657194734889667e-07` pu.  Both equilibrium residuals are
small, but the load and device Jacobian blocks depend on the endpoint used.

The first nonzero block in the frozen bus ordering is the bus-3 diagonal
`2x2` constant-power-load block, with Frobenius norm
`3.6574351814281946e-06`.  At `lambda_c`, its residual is

```text
[[-2.1086214321286434e-06, -1.4975245505866042e-06],
 [-1.4972754343034467e-06,  2.1085858428193660e-06]]
```

The same voltage-endpoint mismatch changes other load blocks (relative
differences are about `4.223e-7`) and unchanged SG linearizations.  The largest
diagonal residual is the common SG block at bus 39,
`7.025098585394785e-05`; it is not a replacement bus.  Thus these omitted
changes are outside the declared four-bus update support.  The core replacement
device changes themselves are correctly represented by the action update.

## Alternatives ruled out

- **Load physics / parameters:** unchanged; both cases use the same power-load
  data and matched-Q policy.
- **Network block:** Ybus is exactly equal between cases.
- **Ordering/permutation:** bus IDs and index maps are exactly equal.
- **Sign convention / normalization:** the raw network Schur and candidate-port
  Schur agree at the already-passing candidate-port tolerance; the observed
  error is localized to omitted diagonal blocks, not a global sign or basis
  transform.
- **Equilibrium/Jacobian provenance:** confirmed as the cause.  The load and
  device blocks are finite-differenced at two separately converged `x,z`
  endpoints rather than one shared matched operating point.

## Minimal permitted correction

Construct only the *base counterfactual port operator used by M1A* at the
flagship's already-solved H4 equilibrium point: use the flagship states for
unchanged, identical SG devices; initialize the replaced base SG devices at the
same bus voltage and the unchanged base dispatch; and evaluate the unchanged
network/load blocks at that same voltage vector.  Keep the base and flagship
canonical DAE solves, equations, parameters, equilibrium outputs, P4, and
`lambda_c` untouched.  This repairs the construction point of the declared
rank-eight identity; it is not a new physical operating point or model change.

An in-memory diagnostic prototype (no source edits) reduced the base+deltaY
relative residual to `4.3277245574459105e-16`, with scalar `|h(lambda_c)|`
`8.123673114402059e-11`, Schur-root distance `3.2014570199586145e-10`, and
complement condition `4.551908942115789`.  Its action-space backward error was
`2.4185882594213892e-11`, above the frozen `1e-12` gate.  Therefore this
prototype is diagnostic evidence only, not a PASS; the complete official
post-correction ladder must determine the ticket outcome without changing any
threshold.

No implementation code had been modified when this root-cause record was
written.

## Pre-correction regression fingerprints

Computed from the unchanged frozen H4/P4 Python case before the implementation
edit:

- model/configuration fingerprint: `5bd6c2c82d876e9667ba177eff952b958902967f19b887fb1e3d7aba06e4fb5c`;
- equilibrium `(x,z)` fingerprint: `18adf2938aa8263bbdb97165dceba885d7323a935664a1d18f16fbac7294ac34`;
- canonical reduced-DAE `A` fingerprint: `da25d45096f80d872ef40caa968fa2c996ef0cc6b1bb92420ad2869c7cc92751`;
- H4: `86` states, `alpha=0.12700646782830968 s^-1`, `f=0.6222796695036534 Hz`;
- frozen F1 ledger: 16 rows total; all 15 non-H4 portfolios (empty/base
  portfolio plus the 14 nonempty proper subsets) are `STABLE`, while H4 is
  `UNSTABLE`;
- existing exact GFL11 Julia parity ledger:
  `PASS_GFL11_EXACT_REPRODUCTION`, 86 states in both codes,
  `|delta alpha|=4.9398755e-10 s^-1`,
  `|delta f|=1.50509e-11 Hz`.

These are read-only regression baselines; the frozen subset lattice and Julia
campaign were not rerun as part of IAS26-010.
