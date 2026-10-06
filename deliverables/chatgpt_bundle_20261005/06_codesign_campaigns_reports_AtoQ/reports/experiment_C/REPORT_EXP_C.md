# Experiment C — Real IEEE-39 Graph Dynamic Self-Energy

## 1. Executive result

**EXP_C_STATUS: PASS**

This run uses the frozen bus-33 Experiment-A `A_reduced`/state map as its primary real input, audits it against a same-model reconstruction (relative discrepancy 0.0), and rebuilds the all-SG baseline plus nominal single-GFL replacements at buses 30, 35, and 37. Controller parameters are unchanged.

## 2. Relation to Experiments A and B

ExpA input status: **PASS**. ExpB-pre status: **PASS**. ExpB's tested generalized graph basis, exact Schur tools, and synthetic self-energy special case were reused; no A/B output was edited.

## 3. Why Sigma(s) is not globally available

The condensed block at zero remains singular. ExpC does not evaluate `Pi_q(0)`, add an epsilon, or use a pseudoinverse. It works directly with the exact generalized object `Psi(s)=sD0+(L0-LG)+Pi_q(s)`.

## 4. Generalized BND operator

The construction is `Tq(s)=s²M+LG+Psi(s)`. Across all five cases, the sample-point median/p95 algebraic reconstruction errors are 3.137097402395099e-17 / 1.5454154575947225e-16. The analytic `Pi_q'` derivative check has median/p95 error 5.588712197100199e-11 / 4.223080096247069e-10; coordinate conversions and units are in TABLE_C20.

## 5. Zero-frequency condensed singularity

Classification: **PARTIALLY_EXPLAINED**. At relative tolerance 1.0e-10, `nullity(A_cc)=1`; the tested `A_cc²` nullity and threshold sensitivity are in `ZERO_FREQUENCY_NULLITIES.csv`. Since nullities do not agree robustly across tolerances, no semisimplicity or residue claim is made. The fitted empirical low-frequency norm slope is 1.000032247027082, over the predeclared conditioning-limited range. Dominant named right-null states: VIndex(33, :gfl₊pll₊Δω_i_rad_s).

## 6. Graph-backbone definition and canonicality

Primary `LG` is the preregistered M-normalized Hermitian part of `L0`, projected onto the relative-angle subspace. The secondary backbone is the Euclidean gauge-preserving Hermitian part of `L0`. Neither was fitted to pole data.

## 7. Generalized graph basis

For C33: minimum eigenvalue of `M`=0.0026525823848649226, M-orthogonality residual=1.7319027339675558e-14, diagonalization residual=1.7378700631082228e-11, zero modes=[2], negative synchronizing-backbone modes=1. These eigenvalues are not called graph frequencies when negative.

## 8. Exact graph-modal transformation

The frequency-grid p95 of `PhiᴴTqPhi - (s²I+Lambda+Psihat)` is 3.589464951382952e-15. `Psi_hat` is reported as `s D0_hat + DeltaL_hat + Pi_hat`; source and diagonal/off-diagonal norms are in TABLE_C06.

## 9. Generalized graph dissipative operator

The reported operator is `D_G_eff=(1/omega) Im_H{Psi_hat(j omega)}`. The direct harmonic power identity maximum error across cases is 5.698518231827367e-16; the synthetic ExpB Sigma-special-case recovery error is 1.0540250939223871e-16. Negative eigenvalues indicate an indefinite operator under this convention, not instability.

## 10. Real IEEE-39 modal mapping

The full finite poles are from the actual PowerDynamics reduced matrices. The analysis set includes the non-gauge spectral-abscissa pair, four low-damping pairs in 0.05–5 Hz, two additional high-retained-participation pairs, and a PLL-state-dominant pair when available. Results use M-weighted right-eigenvector graph projection, not classical participation factors. See TABLE_C07.

## 11. Exact intermodal Schur self-energy Gamma

`Gamma_k=-Tkr*Trr\Trk` is computed with factored solves. The C33 maximum scalar Schur residual is 0.0.

## 12. Exact pathway decomposition

At each successful diagonal root and selected full-system pole offset, the exact pathway matrix `G_k[l,m]=-psi_kl*(Trr\)[l,m]*psi_mk` is saved with its evaluation point and reference ID; its entries sum to Gamma. Top C33 diagonal pathways for spectral-abscissa pole eig92 (s=-0.38197482732380994 + 3.9342252302154828im, offset 0.0001 + 0.0001im): G_2[10,10]:4.363105270555676e6, G_2[1,1]:2.951436182625239e6, G_2[3,3]:14734.411435436232, G_2[4,4]:427.07742234272666, G_2[5,5]:410.0580492606223, G_2[6,6]:101.84090992820916, G_2[8,8]:83.07801036291404, G_2[9,9]:76.16130955553446, G_2[7,7]:10.721758085463437. Individual pathways depend on the chosen basis, especially inside repeated eigenspaces.

## 13. Uncoupled graph-modal roots

Roots of `s²+nu_k+psi_kk(s)=0` were searched with actual pole seeds, analytic derivatives, line search, and a controller-pole conditioning guard. Only converged roots in the configured 0.05–5 Hz band are retained; failed searches are not filled in. See TABLE_C09.

## 14. Pole-shift predictor

For simple nondegenerate roots, `Delta_s=-Gamma_k(s0)/(2s0+psi'_kk(s0))` is compared to actual full PowerDynamics poles. It remains an asymptotic/empirical approximation.

## 15. Predictor validity regime

Classification: **SUPPORTED**; matched mode count=24, median relative error=0.000905638706335024, max=0.13749974048575472. Eta-bin counts/errors are in RESULTS_EXP_C.json and TABLE_C10. The matching threshold and eta bins were fixed in the config before the analysis.

## 16. Real intermodal susceptibility

The approximate ranking uses `|psi_kl psi_lk|/|t_l0|`, compared with the exact diagonal pathway magnitudes. C33 mean Spearman correlations: coupling=0.8534313725490196, susceptibility=0.5871130030959752.

## 17. Complementary-mode resonance / detuning

Gate C8: **NOT_OBSERVED**. No resonance is manufactured. TABLE_C12 reports coupling product, dynamic detuning, approximate susceptibility, and exact diagonal-path size separately.

## 18. Baseline SG versus bus-33 GFL

The two cases use case-specific graph bases because the retained second-order metric changes when the PLL angle-rate coordinate replaces the SG speed coordinate. Modes are aligned by normalized M-reference overlap and physical q-shape overlap; mode indices are not assumed to match. See TABLE_C13.

## 19. Cross-bus 30/33/35/37 analysis

C30, C33, C35, and C37 reuse the frozen nominal SimpleGFLDC parameters. Baseline C0 is included in TABLE_C15.

## 20. Backbone sensitivity

Status: **PARTIALLY_ROBUST**. The primary and secondary choices, principal angles, pole-mode alignment, and Gamma changes are reported in TABLE_C14 and BACKBONE_AUDIT.md.

## 21. Psi component ablation

Operator-only diagnostics compare `sD0`, `sD0+DeltaL`, and full `Psi`; these are not claimed to be physically realizable plants. Results are in TABLE_C17.

## 22. Gate summary

| Gate | Status | Observed |
|---|---|---|

| C0 input integrity | PASS | A=true, B=true, Ared_rel=0.0 |
| C1 generalized Psi exactness | PASS | 3.137097402395099e-17 / 1.5454154575947225e-16; 5.588712197100199e-11 / 4.223080096247069e-10 |
| C2 zero-frequency singularity | PARTIALLY_EXPLAINED | nullity=1, semisimple=false, p=1.000032247027082 |
| C3 graph backbone and basis | PASS | Mmin=0.0026525823848649226, eM=1.7319027339675558e-14, eL=1.7378700631082228e-11 |
| C4 exact graph transform | PASS | p95=3.589464951382952e-15 |
| C5 exact Gamma Schur identity | PASS | p95=0.0 |
| C6 exact pathway reconstruction | PASS | 5.747216055983061e-11 |
| C7 real predictor | SUPPORTED | SUPPORTED, n=24, median=0.000905638706335024 |
| C8 real detuning mechanism | NOT_OBSERVED | NOT_OBSERVED |
| C9 generalized dissipative power identity | PASS | 5.698518231827367e-16; special=1.0540250939223871e-16 |
| C10 backbone robustness | PARTIALLY_ROBUST | PARTIALLY_ROBUST |
| C11 cross-bus reproducibility | PASS | 5/5 |
| C12 ready for ExpD | YES | true |

## 23. Exact results

Exact: the generalized operator identity, the congruence transform, scalar Schur self-energy, exact pathway reconstruction, and the harmonic power identity under the declared convention, subject to their reported numerical gates.

## 24. Approximate results

Approximate: uncoupled diagonal roots, pairwise susceptibility, and first-order pole displacement. The predictor status is `SUPPORTED`.

## 25. Unsupported / blocked claims

No optimal Kp/Ki; no optimal rho; no H4 claim; no global transient-stability theorem; no EMT claim; no field claim; no universal threshold on chi_comm; no proven chi_G stability certificate; no basis-invariant pathway claim.

## 26. Scientific interpretation

The real-case result is that the exact Schur reduction can be expressed through a generalized dynamic self-energy without requiring a finite zero-frequency static correction. Whether this yields useful low-order physical explanation is case- and backbone-dependent; the report retains measured predictor and detuning outcomes rather than transferring synthetic ExpB conclusions.

## 27. Decision for Experiment D

Ready for ExpD: **YES**. Required exact-input, graph, Gamma, pathway, and cross-bus gates are summarized above. This decision does not authorize optimization in ExpC.

