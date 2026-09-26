# PD39 deviation log

## D-000 — pre-preregistration candidate correction

- Date: 2026-09-13
- Stage: model qualification, before campaign preregistration was frozen
- Issue: an initial implementation assumed bus 39 was the slack bus while inspecting the official table.
- Resolution: the official PowerDynamics IEEE-39 CSV was checked directly; bus 31 is `Slack`, bus 39 is `unctrld_machine_load`. The implementation and all preregistration text were corrected to use slack 31 and candidates `{30,32,33,34,35,36,37,38}`.
- Campaign impact: no campaign-level stability, weakness, or portfolio result had been generated under the incorrect set. The qualification output was regenerated after the correction.
