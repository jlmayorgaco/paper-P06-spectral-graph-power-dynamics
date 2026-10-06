# Feedback-cycle falsification campaign, full IEEE-39 SG/GFL model (2026-10-05)

Question: does the controller-network feedback-core mechanism predicted by an external Python reduced model survive in the existing full model, and can sparse PLL retuning break it?
Order of work: `PREREGISTRATION.md` (hashed) -> baseline parity -> all-SG -> continuation -> delay -> return Q -> cycles -> authority -> repair -> forced response -> events.
Folders: `src/` code, `raw/` machine output, `derived/` tables, `figures/`, `logs/`. Historical directories are read-only.
Claims and their labels: `CLAIM_LEDGER.md`. Failed or negative results: `NEGATIVE_RESULTS.md`.
Linear analysis uses the exact-exponential exported linearisation of the PowerDynamics-based ReducedDAE model (not a new model); nonlinear checks use `DelayedEvents.jl`.
