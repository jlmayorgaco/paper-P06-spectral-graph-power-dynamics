# P1 — frozen GFL11 Julia/PowerDynamics parity

status: PASS
evidence_class: CANONICAL_PYTHON_VS_JULIA_DEVICE_PARITY
canonical_source_manifest: raw/gfl11/canonical_source_manifest.json
equation_source: canonical ieee39_devices.py + 20260911_GFL_REPRODUCTION_SPEC.md section 4.2
cases: 64
max_state_relative_error: 5.136675344253906e-16
max_terminal_current_relative_error: 2.7755575615628914e-17
max_system_base_power_relative_error: 0.0
target_median_relative_error: 0.0
target_max_relative_error: 5.136675344253906e-16
injector_interface: PASS
compile_bus_current_source: PASS
compile_message: compile_bus(MTKBus(GFL11Injector); current_source=true)
infinite_bus_harness: PASS
infinite_bus_residual: 4.217554207116368e-9
infinite_bus_state_count: 15
infinite_bus_repeat_delta: 0.0
infinite_bus_spectrum_count: 15
infinite_bus_spectrum_delta: 0.0
infinite_bus_message: GFL current-source -> LoopbackConnection -> dynamic shunt/network bus -> PiLine -> VδConstraint slack
network_sensitivity_csv: raw/gfl11/p1_network_sensitivity.csv
network_residual_preregistered_limit: 1e-6
network_state_gate: >=11
equilibrium_residual: 2.092721098805179e-13
transfer_frequency_range_hz: 0.01–100
transfer_points: 602
transfer_sigma_minimum: 0.04014798437054641
transfer_condition_maximum: 66695.22931098416
transfer_csv: raw/gfl11/p1_transfer_julia.csv
canonical_transfer_csv: raw/gfl11/p1_transfer_canonical_python.csv
transfer_comparison: raw/gfl11/p1_transfer_comparison.json

The transfer is reported over the full requested frequency range. The conditioning diagnostics are evidence, not a pass criterion; frequencies near singularity must be excluded from any claim of relative transfer accuracy.
