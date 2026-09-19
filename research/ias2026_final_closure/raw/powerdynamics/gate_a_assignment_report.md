# PowerDynamics Gate A independent modal assignment

status: PASS
method: scipy.optimize.linear_sum_assignment on pairwise complex-eigenvalue distance
reference_path: nonmutating/default
reference_tolerance: 1.0e-12
lexicographic_entrywise_comparison: not used for this reconciliation
threshold: maximum matched complex-eigenvalue distance <= 1e-6
source: raw/powerdynamics/gate_a_spectra.csv
result: raw/powerdynamics/gate_a_assignment.csv
