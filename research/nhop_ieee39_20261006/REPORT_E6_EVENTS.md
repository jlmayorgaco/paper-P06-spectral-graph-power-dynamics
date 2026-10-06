# E6 five-event nonlinear campaign (frozen graph)

Runner `src/nhop_events.jl` (unchanged), designs unchanged, all 7 requested TOML files existed. Events: (bus, delta) = (8,-100), (16,+100), (16,-100), (29,+100), (29,-100). Pass guard (runner): F <= 0.5 Hz, RoCoF <= 0.5 Hz/s, 0.9 <= V <= 1.1, limiter-fraction ("slack") >= 0.002. Data: `derived/E6_events.csv`, `raw/events/events_<label>.csv`.

| label | tau ms | rho | n | passed | worst F Hz | worst RoCoF Hz/s | min slack | integration errors |
|---|---|---|---|---|---|---|---|---|
| frozen_n0_tau52_rho875 | 52 | 0.875 | 0 | 4/5 | 0.4911 | 0.1677 | 0.00161 | 0 |
| frozen_n1_tau52_rho875 | 52 | 0.875 | 1 | 5/5 | 0.4706 | 0.2072 | 0.00753 | 0 |
| frozen_n3_tau52_rho875 | 52 | 0.875 | 3 | 0/5 | 0.4587* | 0.1775* | 0.0985* | 3 |
| frozen_n0_tau48_rho875 | 48 | 0.875 | 0 | 5/5 | 0.4811 | 0.1846 | 0.00446 | 0 |
| frozen_n1_tau48_rho875 | 48 | 0.875 | 1 | 5/5 | 0.4718 | 0.2017 | 0.00716 | 0 |
| frozen_n0_tau44_rho900 | 44 | 0.900 | 0 | 0/5 | 0.5852 | 0.1933 | -3.1e-7 | 0 |
| frozen_n1_tau44_rho900 | 44 | 0.900 | 1 | 0/5 | 0.5832 | 0.1968 | -3.8e-7 | 0 |

*n=3 extrema only over the 2 events that integrated; the other 3 are not included.

## Every failure
- n0 tau52: event 2 (bus 16, +100): slack 1.61e-3 < 0.002 (F 0.462, other metrics fine).
- n3 tau52: event 1 (bus 8,-100): voltage band, Vmin 0.9282, Vmax 1.1224 (> 1.1). Event 3 (bus 16,-100): Vmin 0.9299, Vmax 1.1215. Events 2, 4, 5: "DDE integration MaxIters" (integration errors, no metrics).
- n0 tau44 rho0.90: all five events F > 0.5 Hz (0.5155 to 0.5852); events 2 and 4 also slack slightly negative (-2.6e-7, -3.1e-7).
- n1 tau44 rho0.90: all five events F > 0.5 Hz (0.5138 to 0.5832); events 2 and 4 slack negative (-3.8e-7, -2.7e-7).

## Does larger n change nonlinear performance?
Not in a helpful way, and not monotonically. n=1 changed little versus n=0: tau 52 improves from 4/5 to 5/5 (only because the limiter slack of event 2 rises from 0.0016 to 0.0075, a marginal case), tau 48 stays 5/5 with ~0.01 Hz lower F and slightly higher RoCoF (0.185 to 0.202), and tau 44 rho 0.90 fails identically (F ~0.58 Hz) for both. n=3 at tau 52 is clearly worse: 0/5, with voltage excursions outside 0.9-1.1 and three DDE integration failures, despite the linear design being stable (E3). Only n=0 vs n=1 was tested for tau 48 and 44; n=2 and n=3 there were not run. Linear stability from E3/E4 does not predict nonlinear success.
