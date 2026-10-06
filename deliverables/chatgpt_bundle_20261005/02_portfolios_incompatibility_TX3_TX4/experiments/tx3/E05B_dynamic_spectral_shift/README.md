# TX3 E05B — Dynamic spectral-shift decomposition

This stage separates the accepted E03 full-descriptor externality into algebraic,
descriptor-mass, and finite-pole factors, then evaluates the finite-pole component on
a new 16-point holdout. It also constructs the connected resolvent trace `Xi_S(s)` and
single-pole contour moments. It does not reopen the negative E05 damping-materiality
decision and does not authorize E05C or E06.

Reproduction order:

1. `preregister_e05b.py`
2. `run_e05b_decomposition.py --population development --workers 8`
3. `freeze_e05b_candidates.py` and commit the freeze
4. `run_e05b_decomposition.py --population holdout --workers 8`
5. `run_e05b_contours.py`
6. `finalize_e05b.py`
7. `validate_e05b.py`
8. `package_e05b.py`
