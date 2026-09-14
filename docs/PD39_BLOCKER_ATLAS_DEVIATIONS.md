# PD39 Blocker Atlas v1 — deviations

This file is frozen before new numerical results. It must be appended only;
existing entries must never be edited or deleted.

| date | item | preregistered rule | deviation | reason | impact |
|---|---|---|---|---|---|
| 2026-09-14 | TDS runtime guard | Same frozen Rodas5P setup to 20 s | Added `maxiters=100000` and explicit retcode recording after a selected C12 blocker exceeded the integrator's default iteration budget; no pulse, duration, tolerances, observables, or acceptance threshold changed. | Prevent unbounded solver time on the preregistered large-alpha stress case. | A max-iteration return is classified as TDS failure, never as physical confirmation; partial traces are retained when available. |
