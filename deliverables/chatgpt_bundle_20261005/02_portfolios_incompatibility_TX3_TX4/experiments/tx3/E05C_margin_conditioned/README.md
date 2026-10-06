# TX3 E05C — Margin-conditioned connected pole externalities

E05C tests whether all 15 E05B-reproducible finite-pole coalitions consume a
material fraction of lower-order modal margin under one independently frozen
uniform constant-power-factor load-stress coordinate. It uses eight new Sobol
seed paths and nine normalized stress levels. E05 and E05B remain immutable.

Execution order:

1. `preregister_e05c.py`
2. commit the preregistration
3. `calibrate_e05c_stress.py` (EMPTY coalition only)
4. `freeze_e05c_endpoints.py` and commit the endpoint freeze
5. `run_e05c_campaign.py --workers 8`
6. `finalize_e05c.py`
7. `validate_e05c.py`
8. `package_e05c.py`

No E06 or new ParaEMT scientific run is permitted by this stage.
