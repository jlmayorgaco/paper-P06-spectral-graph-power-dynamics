# P1 — Fixed-architecture algebraic identities

Status: **PASS**. The evaluated equations are the rank-one PLL pencil update, the rank-two single-port SG retention determinant, and the 4×4 two-port update with cross paths.

- PLL pair: `ΔA = (ΔKp e₄/τ + ΔKi e₅) qᵀ`; E remains identity. All ten GFL sites and two frequencies were checked. Maximum rank residual `1.0082169809793758e-17`, shared-detector relative error `6.572823944327929e-16`, pair affine residual `0.0`, determinant same-device mixed term `2.1317707282671036e-8`, determinant-lemma residual `3.136130363824999e-8`. Near-mode pencil condition reaches `6.007871469016482e11`, which explains the Float64 residual; far-frequency errors are also recorded.
- Exact rational 3×3 check: `(same_device_cross = 0.0, different_device_cross = -30.0, exact = true)`. A deterministic cross-device witness has buses 35,38, `|cross|=0.0012453493326376323`.
- SG: `T=TF+εΠ(YF−YSG)Πᵀ`; single-port operator/determinant errors `7.805322367159452e-18` / `5.106253525807105e-14`. Two-port 4×4 errors `3.871016991356509e-18` / `5.062616992290714e-14`; cross term `0.21236729879516358`.
- PI geometry: ω=0 is handled as a real line (`true`); positive-frequency cases use the real 2×2 solve only when conditioned. No inverse of the singular real-frequency matrix is used.
- Certificate scope: conditional exact identities on the frozen ExpN architecture and tested frequency points. These identities do not establish a global co-design optimum.
- Time: `56.435` s. Next allowed stage: P2 if all measured identity gates pass. Unresolved work: exhaustive PI geometry and continuous joint optimum.
