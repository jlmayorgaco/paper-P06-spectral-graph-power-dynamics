# ParaEMT EMT campaign — preregistration amendments

Each amendment is committed before the runs it affects. Preregistration v1:
commit c2947bd8.

## A1 — zero-sequence reference of the unit-test infinite bus (implementation; affects EMT02/EMT03 only)

**Found.** The first EMT02 run was non-finite within 1 ms, before any
comparison was computed.

**Cause.** The synchronous-frame algebraic Norton elements (prereg §3) have
no zero-sequence path by construction:

`Σ_m a^{m−k} = 0`, so the all-ones vector lies in the null space of their
3×3 blocks.

The preregistered unit-test system has no element to ground at all:
- device at T;
- an R–L line with no charging;
- a Thevenin infinite bus.

The zero-sequence network is therefore floating, and the network matrix is
singular (condition number 3e16). The IEEE-39 cases are not affected, because
the line-charging capacitors and the bus-4/5 shunts provide the zero-sequence
reference to ground (EMT01 equilibrium passed).

**Amendment.** The infinite bus is realized as a **solidly grounded**
Thevenin source. Its neutral grounding is a zero-sequence-only conductance
`G_z = (g_z/3)·J`, where J is the 3×3 all-ones matrix and g_z = |1/z_s| = 1e4.
It is stamped at bus INF.

- Since J·v = 0 for positive- and negative-sequence sets, G_z adds nothing to
  the positive-sequence network.
- The preregistered positive-sequence test system (z_s, the line, the device)
  is unchanged.
- No tolerance, disturbance or pass rule changes.
