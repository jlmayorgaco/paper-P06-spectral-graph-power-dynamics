# 1. Executive summary (placeholder — finalize after all phases complete)

This section is filled in at the end (Phase 25), after every gate is known.
Below is a running draft, updated as phases complete; do not treat as final
until the campaign's last commit.

## Confirmed so far (E1, E1b, E10, E3, E4)

- **GOLD-A components (E1, E2 pending):** A1 and A3 PASS on 18/18 base-stable
  holdout policies — contextual sign reversal recurs, and even the
  best-possible fixed node ranking is materially insufficient (median regret
  rate 37%). A2 (mode-tracked) also PASSES.
- **HS:** α_⊥ is neither submodular nor supermodular anywhere tested (39/39
  policies), and the tracked-mode margin is not more structured. NEGATIVE BUT
  INFORMATIVE.
- **GOLD-B: PASS.** Total dynamic line sensitivity (Dtotal) beats the best
  static baseline by 0.49 Spearman (0.987 vs 0.50), with 0.90 top-5 precision.
  This is the strongest result so far.
- **H2 (frozen ≠ total):** TRUE for both nodes and links, under both
  equilibrium semantics. The implementation-validity gate is essentially
  perfect (8002/8002).
- Shapley/context-averaged diagnostics are **not useful** on this class: they
  are dominated by rare, physically real fast-instability contexts.

## Pending (E2, E6, E7, E8, E9, E11-E14, E16, E23; E17 Africano gate)

Updated below as each completes.
