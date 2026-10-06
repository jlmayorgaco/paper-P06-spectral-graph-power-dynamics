# D3 lossless network Taylor validation

Status: NUMERICALLY_VALIDATED for the declared lossless sinusoidal branch-power map only. It does not validate lossy AC Taylor tensors, a nonlinear stability proof, or DDE frontier predictions.

Frozen protocol: 12 seeded zero-mean angle directions, amplitudes 0.0001, 0.001, 0.01, 0.05, 0.1 rad; 60 comparisons total. See TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv.

- first_order_relative_error: max=0.0011084378059901478, median=3.656368598546662e-5.
- second_order_relative_error: max=0.00029633477013230855, median=2.280844928605951e-6.
- third_order_relative_error: max=3.057305316509495e-7, median=5.0512522681735905e-11.
