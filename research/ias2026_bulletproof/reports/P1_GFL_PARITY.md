# P1 — frozen GFL11 Julia/PowerDynamics parity

status: STOPPED_BY_GATE
evidence_class: FRESH_DEVICE_LEVEL_REPRODUCTION
equation_source: 20260911_GFL_REPRODUCTION_SPEC.md section 4.2
cases: 64
max_state_relative_error: 5.136675344253906e-16
max_terminal_current_relative_error: 2.7755575615628914e-17
max_system_base_power_relative_error: 0.0
target_median_relative_error: 0.0
target_max_relative_error: 5.136675344253906e-16
injector_interface: PASS
compile_bus_current_source: PASS
compile_message: compile_bus(MTKBus(GFL11Injector); current_source=true)
infinite_bus_harness: FAIL
infinite_bus_residual: Inf
infinite_bus_state_count: 0
infinite_bus_message: ArgumentError: reducing over an empty collection is not allowed; consider supplying `init` to the reducer
equilibrium_residual: 2.092721098805179e-13
transfer_frequency_range_hz: 0.01–100
transfer_points: 81
transfer_sigma_minimum: 0.04014798437054641
transfer_condition_maximum: 66695.22931098416
transfer_csv: raw/gfl11/p1_terminal_transfer.csv

The transfer is reported over the full requested frequency range. The conditioning diagnostics are evidence, not a pass criterion; frequencies near singularity must be excluded from any claim of relative transfer accuracy.
