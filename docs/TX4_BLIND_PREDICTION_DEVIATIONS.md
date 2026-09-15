# TX4 blind-prediction deviations

Initialized before blind prediction. No deviations are authorized. Any
deviation must state date/time, frozen rule affected, reason, exact action,
stop/continue decision, and claim impact. Negative outcomes remain in the
final report.

## Recorded deviations

1. **Reduced g-boundary timing (post-reveal, 2026-09-14).** The first blind
   V4/V9 run and its hashes were committed before the full-order reveal. The
   clean g-boundary search was then implemented and executed as a reduced-only
   reproduction, but after the reveal comparison process had already loaded
   historical answer tables. The boundary solver itself does not import an
   answer key and records `protocol_status=REDUCED_REPRODUCTION_AFTER_REVEAL`.
   The numerical boundary agreement is therefore reported as a reproduction,
   not as a fully blind C05 prediction. This downgrades the strict predictive
   claim but does not change the V4 blind portfolio result.

2. **Contextual-return implementation provenance.** The contextual-return run
   reused the audited TX4 closure implementation from the immediately prior
   TX4 worktree and was then materialized on this branch. Its exact source is
   now included as `reports/poster/ias2026/research/experiments/tx4_contextual_return.py`.

3. **V9 speed.** The reduced closure was slower than the direct full-order
   benchmark on this hardware. This is retained as a negative runtime result;
   no speedup claim is made.

Initial status: `RECORDED; NO POST-HOC PARAMETER EXPANSION`.
