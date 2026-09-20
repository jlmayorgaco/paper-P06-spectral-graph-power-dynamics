# P3/L6 — same-model collective mechanism

status: `JULIA_COLLECTIVE_MECHANISM_VALIDATED`
evidence_class: `TRUE_FROZEN_CUSTOM_MODEL_JULIA_PORT_OPERATOR`

The Julia port operator was assembled from the same frozen SG/AVR/PSS/GFL DAE used by the cross-code census. The four-target V4 update is represented in an 8-dimensional rectangular terminal action space.

At the minimum over the frozen 0.3–1.5 Hz imaginary-axis grid, frequency=0.575 Hz, the closest collective Q eigenvalue is -1.0829168901378603 + 0.029752479141036896im, and its distance to -1 is 0.0880932499410246. The minimum local factor singular value over the grid is 0.29424879369829876; therefore the collective near-closure is not explained by a singular local factor.

The determinant identity det(I+M)=det(I+M_local)det(I+Q) is checked at every grid point; the maximum normalized residual is 7.216449660063518e-16. Linear solves use T0\U and the maximum normalized backward residual is 1.9092140494856664e-17.

The CSV records full, individual, collective, Q, local-regularity, and solve diagnostics. The sign-dual diagnostic -Q is retained through the recorded q_closest value; no universal sign claim is made outside this frozen port convention.
