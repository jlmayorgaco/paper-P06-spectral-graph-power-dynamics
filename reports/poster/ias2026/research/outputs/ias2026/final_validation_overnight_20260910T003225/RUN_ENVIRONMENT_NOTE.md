# Run environment note

Recorded during the run, not after.

The machine carried an **unrelated concurrent workload** for part of this run: a
separate `localizer_detector_v8` job with roughly twelve worker processes under a
different interpreter. Worker counts here were chosen from `os.cpu_count()` (22)
and free RAM, and did not account for a foreign workload, so the later
experiments ran oversubscribed — roughly 28 busy processes on 22 cores.

This affects **wall-clock time only**. Every experiment is seeded and
deterministic, all parallelism is over independent samples with no shared state,
and BLAS threads are pinned to one inside every worker. Timings in the manifests
of the later experiments are therefore not comparable with the early ones
(E33 ran in 20 s on a quiet machine; E36 took 792 s for a comparable amount of
work on a busy one).

Nothing was retried, re-seeded or re-tuned because of this.
