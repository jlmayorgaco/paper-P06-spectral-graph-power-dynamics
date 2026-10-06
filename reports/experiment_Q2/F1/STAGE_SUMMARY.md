# F1 — Small-disturbance grid-frequency identity

Status: **PASS_LOCAL_LINEAR_IDENTITY**. Eight configurations (all-SG, corrected ExpG, ExpN, and five seeded mixed designs) were independently trimmed and simulated in PowerDynamics at 1 MW and 5 MW; the analytical reduced model was evaluated at the same points. Passive bus-frequency traces use the fixed 0.1 s SG derivative window and are compared at buses 30–39 from 1.2 to 6 s.

- Maximum normalized linear/PD discrepancy: 0.0010049254311783973
- Maximum normalized 1-to-5 MW scaling discrepancy: 0.0008034184859158516
- Maximum P/Q trim error: 5.732658792112489e-13 pu
- Acceptance gate: 3% for both trace identity and amplitude scaling.
- Result table: `TABLE_Q2_F1_grid_frequency_identity.csv`.
