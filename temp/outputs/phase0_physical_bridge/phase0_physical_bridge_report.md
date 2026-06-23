# Phase-0 Physical ANDES-to-(L,M,D) Bridge Audit

Gate status: **BLOCKED**

Acceptance rule: critical mode frequency error <= 5%, first-six mean frequency error <= 5%, and positive damping sign reproduced by physically extracted scalar D.

## base
- Accepted: False
- Full modes in 0.1--3 Hz band: 11
- Surrogate modes: 9
- Blocking reason: Physical L and M place several frequencies in the electromechanical band, but GENROU.D is zero in the packaged case, so the extracted scalar D cannot reproduce the positive damping of the full ANDES modes. Controller/exciter/PSS states provide damping that is not captured by a scalar nodal D without a validated reduction.
- Mode match CSV: `outputs\phase0_physical_bridge\base_mode_match.csv`
- Matrices JSON: `outputs\phase0_physical_bridge\base_matrices.json`

## no_pss
- Accepted: False
- Full modes in 0.1--3 Hz band: 11
- Surrogate modes: 9
- Blocking reason: Physical L and M place several frequencies in the electromechanical band, but GENROU.D is zero in the packaged case, so the extracted scalar D cannot reproduce the positive damping of the full ANDES modes. Controller/exciter/PSS states provide damping that is not captured by a scalar nodal D without a validated reduction.
- Mode match CSV: `outputs\phase0_physical_bridge\no_pss_mode_match.csv`
- Matrices JSON: `outputs\phase0_physical_bridge\no_pss_matrices.json`

Dependent full-ANDES estimator validation phases must remain blocked. The result is a useful negative audit: the scalar network surrogate needs a validated controller-state reduction before it can be used as full-ANDES ground-truth evidence.
