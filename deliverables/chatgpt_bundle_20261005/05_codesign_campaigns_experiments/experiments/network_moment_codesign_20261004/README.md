# Network frequency moments and collective PLL latency compensation

Start with `REPORT_ES.md` for the executed result and `POSTER_CLAIMS.md` for
claim limits. `THEORY_MOMENT_SECTION.tex` contains the full-network proof and
conditional graph-control corollary; it is integrated into the user's existing
`experiments/theory_collective_damping_20261003/THEORY.tex`.

The exact scalar coefficient law connects the synchronization residue of the
physical network to allowable finite PLL gain/delay changes. It distinguishes
low-frequency response matching from finite-frequency damping. It does not
compute maximum GFL replacement. The fixed prior87.6462215% design is an input.

Evidence map:

- TABLE01--04: exported input parity, signed network weights, coefficient
  identities and finite-frequency asymptotic error.
- TABLE05: exact-delay pole derivatives versus centered differences.
- TABLE06--09: preserved initial multimode attempt, including its failure.
- `globalized/TABLE13--15`: final exploratory controller attempt, ablations,
  complete gain vectors and eleven-mode catalog.
- TABLE10--12: independent Julia parity, complete-region DDE contour and
  five nonlinear delayed events for the headline candidate.
- TABLE16 and MOMENT_ASSUMPTION_CERTIFICATE.json: Arb enclosures of hidden
  regularity, simple rotational derivative and network weight signs.
- REVIEW_MATH/EVIDENCE/NOVELTY.md: independent reviews and bounded literature map.
- Protocols, addenda, source locks and failure logs: preserved research history.

Read `REPRODUCE.md` before executing code. `STATUS.json` is the generated
machine-readable result, including solver nonconvergence and publication limits.
The complete ZIP includes figures, trajectories, reports, theory, code, required
read-only repository inputs and an exact-byte inventory. No commit is needed
to reproduce an extracted bundle, but Python/Julia packages must be installed.
