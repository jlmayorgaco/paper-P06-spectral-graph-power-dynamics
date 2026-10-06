# ExpH post-freeze PowerDynamics validation

Candidate SHA-256: `2967b772d0330e71e528298940e02812b0642c0a71610f2a72bd1c16eedaeac3`. The SHA gate passed before PowerDynamics was loaded. The frozen candidate was not changed.

Overall post-freeze status: **PARTIAL**. Candidate equilibrium residual: 1.3951473992034526e-12; PD spectral abscissa: -0.0013257499788513162 s⁻¹; analytical abscissa: -0.06479596417248183 s⁻¹; difference: 0.06347021419363051 s⁻¹.

Candidate has 102 finite raw poles and 3 numerical gauge poles. The comparison is between different model realizations; it is a margin-level check, not a pole-by-pole identity.

Nonlinear load pulses at bus 16: 0.0005 and 0.001. TDS scaling gate: PASS. See TABLE_H21_tds_validation.csv and TABLE_H21_tds_trajectories.csv.

Controller compatibility error for the frozen Kp/Ki under PD39's shared PLL scale is 5.741790543023484e-8.

This validation does not certify the analytical full-block uncertainty bound in the nonlinear PowerDynamics model, nor local/global optimality of the design.
