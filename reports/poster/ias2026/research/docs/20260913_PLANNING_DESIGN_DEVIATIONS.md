# Deviation log: planning and design campaign (PD1, PD2)

Preregistration: `docs/20260913_PLANNING_DESIGN_PREREG.md`. Entries are added as
they occur; no rule is changed after results are visible.

1. (before freezing) PD2 created its checkpoint directory only in the runner, so
   a direct call of `run_task` failed on the partial checkpoint. Fixed by creating
   the directory at task start. Found in a smoke test on a non-campaign input.
