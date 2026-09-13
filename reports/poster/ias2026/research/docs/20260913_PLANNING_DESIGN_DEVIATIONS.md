# Deviation log: planning and design campaign (PD1, PD2)

Preregistration: `docs/20260913_PLANNING_DESIGN_PREREG.md`. Entries are added as
they occur; no rule is changed after results are visible.

1. (before freezing) PD2 created its checkpoint directory only in the runner, so
   a direct call of `run_task` failed on the partial checkpoint. Fixed by creating
   the directory at task start. Found in a smoke test on a non-campaign input.
2. (PD1 execution) The PD1 pool completed 239 of 240 samples in about 23 min;
   the worker that held sample 1596 terminated without a Python exception
   (the pool spawned a replacement worker and `imap_unordered` then waited
   indefinitely for the lost task). After confirming that no PD1 worker was
   computing, the PD1 process tree was stopped (only this campaign's processes;
   a foreign `run68.py` job was left untouched). Sample 1596 was then run alone
   with the same code and `-X faulthandler`; it completed normally in 143 s
   (84 SLSQP iterations, 745 evaluations). Its record is appended to
   `PD1_records.jsonl` with `rerun_isolated: true` and without the solution
   vector `x` (not printed by the isolated run). No rule, target or threshold
   was changed. Log: `results/PLANNING_DESIGN/PD1/PD1_sample1596_isolated.log`.
