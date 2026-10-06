# IAS26-030 — F2 closure and boundary audit

Status: **PASS**

Run ID: `20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1`

F2 was assembled from the frozen baseline only. No equilibrium, eigenvalue, continuation, or root solve was run. The 17 alpha/collective sweep points were exact-joined to their physical pre-normalized local factors in fixed port order `(30, 33, 35, 37)`. The exact stored boundary row was added from the boundary artifact, not interpolated.

- `g*`: `0.20768140519037842`; eigen/return boundary difference: `3.72287e-10` (baseline tolerance `5e-05`).
- At the stored boundary, minimum physical local `sigma_min(I+M_ii)`: `0.489235486`; collective `sigma_min(I+Q_H)`: `4.39367037e-08`.
- Maximum local-versus-collective frequency difference in matched records: `1.63328e-11 Hz` (reported, no new threshold).
- Maximum boundary Schur identity residual: `2.29593e-16`; frozen port identity maximum: `1.43049e-16`.
- Mode audit: alpha monotone=True; frequency monotone=True; no eigenvectors/gaps were retained, so mode-family identity is not certified.

The figure supports only the narrow claim that physical local factors remain regular while collective closure approaches singularity. It does not establish positive local damping or validate reduced poles against the DAE.
