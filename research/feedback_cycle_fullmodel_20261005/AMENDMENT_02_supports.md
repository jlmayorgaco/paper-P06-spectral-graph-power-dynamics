# AMENDMENT 02 (2026-10-05, after Phase H/F, before any Phase J run)
S0 = {30,33,36,37} was tested first and failed in all three variants (see NEGATIVE_RESULTS). The comparator supports below are fixed here by rule BEFORE they are run:
- sens4: the four sites with the largest joint authority sqrt(dRe/dlogKp^2 + dRe/dlogKi^2) at the leading 44 ms root (TABLE08): 30, 37, 39, 38.
- phys4: endpoints of the two physically strongest pairs (TABLE_G01): 33-34 and 35-36 -> {33,34,35,36}.
- core4: endpoints of the two pairs with largest |log|1-p|| at the leading 44 ms root (TABLE_F02): 35-36 and 30-37 -> {30,35,36,37}.
- uniform: all ten sites scaled by one common factor (separate for Kp and Ki); all10: unrestricted nodal gains.
Same engine (src/repair.py), bounds (0.25-4x nominal), sigma_req = 0.05, tau = 44 ms, step cap 0.1, <= 60 iterations, success = full-spectrum count N_margin = 0.
