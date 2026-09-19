# Open Limitations

The following items remain open and are not represented as completed results:

- same-model PowerDynamics implementation of the frozen L0 GFL dynamics;
- the P1 minimal current-source/infinite-bus topology currently retains zero
  network states despite passing device-level Julia/Python parity;
- P2–P6 remain dependency-stopped until that P1 integration issue is resolved;
- a second independent converter model and a genuinely new blind holdout;
- a preregistered physical common uncertainty envelope and robust radius;
- a new Julia nonlinear TDS run against the frozen converter model;
- current-limit, DC-link, EMT switching, protection, hardware, and full repair-path tests;
- a new asymptotic scaling sweep beyond the historical V9 reveal.

These limitations constrain the safe poster and paper claim. They do not erase
the frozen IEEE-39 evidence or the negative V9 transfer result.
