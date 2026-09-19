# Final Bulletproof Audit

## Executive disposition

The corrective campaign closes the claims supported by the frozen evidence and
explicitly refuses to promote the rest. Gate 0 and PowerDynamics Gate A pass
within their declared scopes, but the campaign is not a final paper/poster
release because the dependent same-model, robustness, second-model, holdout,
and new-TDS gates remain incomplete.

## What survived

1. The frozen IEEE-39 P4 nominal H4 alpha reproduces as
   `0.1270064666836382 s^-1` against `0.1270065 s^-1`.
2. The frozen P4 structure reports 15 stable proper subsets, supporting the
   benchmark-specific minimality statement.
3. The randomized theorem/property tests pass across four modules, including
   the port bridge and transverse construction. The contextual-return test
   states and enforces its `Q_ii=0` assumption.
4. The frozen F7A boundary witness and sign audit reproduce retrospectively.
5. The frozen G2 nonlinear table reports 32/32 declared TDS agreements.
6. The documented governor configuration changes P4 H4 alpha from positive to
   negative (`-0.07454098428813066 s^-1`). This is a policy-conditioned
   stabilization/repair observation, not a universal repair certificate.
7. The official PowerDynamics 5.0.0 IEEE-39 tutorial passes Gate A with three
   tolerances, three deterministic initial guesses, two initialization paths,
   residual checks, eigenvalue consistency, and Jacobian conditioning. This is
   package-level evidence only.

## What was weakened or refuted

- The corrected same-policy global V9 scope contains 402/512 correct, 110
  false-safe, and 0 false-unstable cases; all false-safe cases are
  `APERIODIC_REAL`.
- The corrected targeted 0.3–1.5 Hz scope contains 511/512 correct, 1
  false-safe, and 0 false-unstable cases. It is not a global certificate.
- The corrected target-family blocker antichain is exact: 14 blockers, orders
  5x4, 6x5, 3x6, and kappa=4.
- The P4 blocker does not survive the documented governor policy; the sign
  changes in the frozen evidence.
- No universal IEEE-39, all-policy, all-model, or all-uncertainty claim is
  justified.

## What was not established

- PowerDynamics same-model GFL parity or independent reproduction of the
  spectral/contextual-return mechanism.
- A second converter model, a new blind topology/model holdout, a physical
  common uncertainty envelope, a robust LFT radius, or a new Julia TDS run.
- EMT, current-limit, DC-link, protection, switching-path, or hardware
  certification.

## Direct answers to the closure questions

| Question | Closure answer |
|---|---|
| Strongest survived result | Benchmark- and policy-conditioned P4 minimal incompatibility with exact frozen nominal value, theorem identities, and retrospective nonlinear agreement. |
| Strongest weakened result | Transfer/generalization beyond the frozen case; V9 is negative evidence. |
| PowerDynamics role | Gate A validates the official IEEE-39 tutorial only; no same-mechanism claim. |
| Second model | Not executed; no `SECOND_MODEL_VALIDATED` label. |
| Robust radius | Not computed; no physical robustness claim. |
| Collective closure vs. robustness | Frozen nominal closure survives; common-envelope robustness remains untested. |
| Repair | Documented governor reverses the P4 H4 sign; full-path repair is not certified. |
| Unseen case | No new blind unseen-case result; historical V9 is retrospective. |
| Safe poster wording | Use the scoped claim in `docs/CLAIM_LEDGER.md`. |
| Supplementary material | Include tables, logs, frozen-source paths, labels, and all negative/stopped gates. |

## Final safe claim

In the frozen IEEE-39 benchmark and stated policy, the spectral/collective
closure audit reproduces a P4 minimal incompatibility and the declared frozen
nonlinear verdicts. The corrected same-policy V9 transfer is scope-dependent
(402/512 globally; 511/512 in the declared 0.3–1.5 Hz band), the documented
governor changes the P4 H4 sign, and robustness, second-model, same-mechanism
parity, and new blind-holdout claims remain untested or stopped by scope.
