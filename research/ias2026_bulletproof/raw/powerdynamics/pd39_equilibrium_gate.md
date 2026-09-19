# PowerDynamics IEEE-39 Gate A

status: POWERDYNAMICS_VALIDATED
gate: A
scope: official PowerDynamics tutorial model; not SG->GFL parity
package_version: 5.0.0
julia_version: 1.11.9
buses: 39
branches: 46
nonmutating_paths: 9
tolerances: 1.0e-8, 1.0e-10, 1.0e-12
initial_guess_paths: default, guess, perturbed
residual_limit: 1.0e-8
spectrum_consistency_limit: 1.0e-6
mutating_componentwise_status: PASS
mutating_componentwise_residual: 7.966960424710123e-13
mutating_componentwise_spectrum_delta: 7.138278946393519e-9
paths_csv: raw/powerdynamics/gate_a_paths.csv
source: PowerDynamics official docs/examples/ieee39_part1.jl
