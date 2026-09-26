# Negative Results Ledger

This file is append-only for the new campaign. A STOP is a scientific result.

| date | phase/gate | result | consequence | evidence |
|---|---|---|---|---|
| 2026-09-18 | environment | Julia 1.11.9 is present; `PowerDynamics` was absent globally but installed in the isolated campaign project | The official tutorial equilibrium gate passes; same-model GFL parity remains `STOPPED_BY_GATE` | `raw/powerdynamics/pd39_equilibrium_gate.md`, `docs/POWERDYNAMICS_RECONCILIATION.md` |
| 2026-09-18 | environment | Runtime Python is 3.13.14 with NumPy 2.3.5/SciPy 1.16.3, while the frozen record reports NumPy 2.5.2/SciPy 1.18.1 | Record as environment deviation; do not silently call this exact environment reproduction | `reports/REPOSITORY_INVENTORY.md`, `docs/DEVIATIONS.md` |

Historical negative evidence is preserved and linked in `docs/REVIEWER_ATTACK_MATRIX.md`; it is not deleted or rewritten here.
