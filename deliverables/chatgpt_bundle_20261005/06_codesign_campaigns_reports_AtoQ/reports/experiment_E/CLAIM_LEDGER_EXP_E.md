# ExpE claim ledger — current evidence

| Claim | Status | Evidence and limit |
| --- | --- | --- |
| Ten generators, buses 30–39, are in the frozen independent design domain | Exact data discovery | E01, E02, frozen TOML and SHA. |
| Initialized SG dispatch is 5402.761 MW | Analytical from frozen machine states | E01, E03; nominal load setpoints give a different and inapplicable total. |
| The analytical all-SG model reproduces ExpC | Numerical identity within roundoff | Spectral Hausdorff 2.05e-12. |
| One-bus replacements reproduce archived ExpC cases | Numerical identity within roundoff | E04; maximum spectral distance 1.76e-11. |
| Independent PLL gains admit a one-scalar Woodbury port | Exact algebra and numerical check | E05; median/p95 relative error 2.64e-16/4.81e-16. |
| Multi-bus port closure has 20 coordinates and satisfies determinant lemma | Exact algebra and numerical check | E06; maximum sampled logdet error 1.43e-13. |
| All-GFL endpoint has an additional near-zero frequency mode beyond the angle gauge | Strong numerical evidence; formal all-gain proof open | E07 at nominal and two independent gain patterns. |
| 5401.571878 MW can be reached in the analytical model at the required margin | Verified analytical provisional point | E13 step 7, E14, E15; full spectrum and KCL checked. |
| This point is the maximum over the frozen 30-variable domain | **Unverified** | Branch enumeration, active-set KKT/SOSC, global gap remain open. |
| The provisional point agrees with independent detailed-model eigenvalues and TDS | **Unverified** | No ExpE validation run; final candidate has not been frozen. |
