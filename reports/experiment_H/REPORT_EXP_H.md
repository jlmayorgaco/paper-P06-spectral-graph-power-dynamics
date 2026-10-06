# Experiment H — Jordan / self-energy SG→GFL robust co-design

## Executive result

**Overall theory result: PARTIALLY_SUPPORTED.** The all-GFL quotient structure, linear retained-SG root scaling, gain-dependent authority, exact Woodbury port form, collective self-energy cancellation, graph-subspace alignment, and analytical direct small-gain certificate are supported by the frozen analytical model. The final full-closure optimum is not certified, PLL gains remain essentially nominal, the transient frequency diagnostic is active and was not fed back into the KKT problem, and the independent PowerDynamics model fails the requested spectral margin and disagrees materially with the analytical abscissa.

The frozen analytical robust candidate retains 11.2 MW of SG at bus 36 (\(\epsilon=0.02,\rho=0.98\)); it converts 5391.561 MW (99.7927%) to GFL. Its analytical spectral abscissa is \(-0.064796\,s^{-1}\), below the nominal \(-0.05\,s^{-1}\) boundary. At \(\beta=10^{-4}\), the direct resolvent peak is 7748.518 at 0.04510 rad/s, giving small-gain margin 0.22515. Candidate SHA-256: `2967b772d0330e71e528298940e02812b0642c0a71610f2a72bd1c16eedaeac3`.

## Isolation and reproducibility

Work was performed on branch `research/expH-jordan-mixed-codesign`. The preblind predictor was frozen and hashed before reading ExpE/ExpG final candidates. The final candidate was frozen before loading PowerDynamics. The analytical design entry point and `src/bnd_design_h/` do not import PowerDynamics. No candidate field was changed during post-freeze validation. No push was made. Artifacts A–G were not edited as part of ExpH.

The frozen domain discovers buses 30–39 programmatically and reproduces initialized SG dispatch of 5402.761089978776 MW. The analytical network preserves initialized ZIP loads, including the bus 31 and bus 39 data. `tables/TABLE_H00_model_integrity.csv` passes archived spectrum and closure checks.

## Hypothesis results

| Hypothesis | Result | Evidence |
|---|---|---|
| H1 — physical zero after gauge quotient | **SUPPORTED** | Nullity 1 / algebraic multiplicity 1 after the explicit gauge quotient; full chain length 2. |
| H2 — gains alone remove the endpoint zero | **FALSIFIED** | The physical quotient zero remains at zero for nominal, three deterministic gain patterns, and both gain-box corners. |
| H3 — retained-SG authority depends on gains | **SUPPORTED** | Mixed derivative validation median (1.28\times10^{-6}), p95 (1.78\times10^{-6}); normalized Jacobian rank 1. |
| H4 — local predictor improves exact boundary | **PARTIAL** | Linear exponent is verified, but exact correction changes retention materially and rejects the predicted best anchor. |
| H5 — Woodbury reduces gain structure | **SUPPORTED** | Exact one-scalar PLL port factorization per device; active mixed-authority subspace is rank 1. |
| H6 — collective self-energy matters | **SUPPORTED** | Direct and self-energy terms strongly cancel; net coefficient would be misrepresented by the direct term alone. |
| H7 — noncommutative graph structure explains active gain direction | **SUPPORTED** | Degree-two graph words capture 99.24% of active gain-subspace energy; principal angle 5.00°. |
| H8 — direct small-gain robustness improves on ExpG | **SUPPORTED analytically** | 5391.561 MW converted at \(\beta=10^{-4}\), versus ExpG's 1580.990 MW at \(\beta_{design}=1.699\times10^{-6}\). Independent PD margin still fails. |
| H9 — true \((\rho,K_p,K_i)\) co-design | **PARTIAL** | A mixed support/controller search was run, but the final gains are effectively nominal and no joint KKT certificate exists. |
| H10 — robust converted MW exceeds ExpG | **SUPPORTED analytically** | +3810.571 MW under the analytical direct-resolvent model. This does not transfer to the independently linearized PD model. |

## Root law and mixed derivative

The explicit quotient removes the rotational direction before classifying the physical root. The full endpoint has nullity one and algebraic multiplicity two; after quotienting, the physical root is simple. Factoring the gauge root from the exact closure determinant gives

\[
F_i(s,\epsilon_i,K)=a_2(K)s+a_{1,i}(K)\epsilon_i+O(s^2,s\epsilon_i,\epsilon_i^2),
\quad \lambda_{phys,i}=-A_i(K)\epsilon_i+O(\epsilon_i^2),
\quad A_i=\operatorname{Re}(a_{1,i}/a_2).
\]

The full-window exact root exponent is 1.01088 on average; all buses are within the predeclared 5% tolerance. At \(\epsilon=0\), gain changes do not remove the zero, while the coefficient (A_i(K)) has nonzero mixed derivatives. Analytic determinant-trace derivatives pass independent complex-step validation. Near-zero derivatives are excluded using the recorded fixed cutoff.

Details are in `DERIVATION_H_JORDAN_ROOT_LAW.md`, `DERIVATION_H_WOODBURY_AUTHORITY.md`, and tables H01–H09.

## Predictor, exact correction, and robustness

The local predictor ranks bus 38 first, at 7.380 MW retained. The exact nominal all-pole correction shows a 78% retention correction on the predicted bus-38 row, and that candidate fails another pole. The best robust candidate found on the fixed single-anchor epsilon grid uses bus 36 at \(\epsilon=0.02\), retaining 11.2 MW. This satisfies the exact analytical full-spectrum margin \(-0.064796\,s^{-1}\) and the declared direct small-gain condition at \(\beta=10^{-4}\). The frontier's estimated \(\beta_\star\) for this candidate is \(1.29057\times10^{-4}\).

The search is deterministic and enumerates ten single-anchor buses at fixed 0.02 epsilon increments with exact spectrum and resolvent checks. It is not a continuous joint optimizer, does not exhaust all pair supports, and supplies no KKT, LICQ, SOSC, local, or global optimum certificate. Do not interpret the candidate as the maximum feasible replacement over the full admissible set.

## Transient screen and independent validation

For the analytical final candidate, the 100 MW transient screen predicts 0.4177 Hz/s peak RoCoF against a 0.5 Hz/s limit, but 3.327 Hz frequency excursion against a 0.5 Hz limit. The frequency constraint is active (`TABLE_H18_transient_active_constraint_test.csv`) and would require reopening the design problem; that KKT correction was not completed. This is a design limitation, not a pass.

After the final-candidate SHA check, the separate PowerDynamics validator reached an equilibrium with residual infinity norm \(1.40\times10^{-12}\), 102 finite raw poles, and 3 numerical gauge poles. It found \(\alpha_{PD}=-0.001326\,s^{-1}\), which does not meet the \(-0.05\,s^{-1}\) requirement and differs from the analytical abscissa by 0.06347 s⁻¹. The fractional bus-36 replacement and candidate PLL settings were representable in the PD39 weighted model; the shared-scale gain mismatch was (5.74\times10^{-8}). Thus this is a measured model disagreement, not a blocked fractional-capacity implementation.

The nonlinear load-pulse TDS at bus 16 passed the 2:1 response scaling screen. Pulses 0.0005 and 0.001 produced peak retained-SG COI deviations 0.0001931 and 0.0003861 Hz, peak RoCoF 0.002892 and 0.005783 Hz/s, and maximum voltage deviations (2.50\times10^{-5}) and (4.99\times10^{-5}) pu. The min/max bus voltage was 0.98196–1.06351 pu. TDS scaling passing does not repair the spectral-margin disagreement.

The separate post-freeze results are in `POSTFREEZE_VALIDATION_EXP_H.md`, `tables/TABLE_H20_powerdynamics_validation.csv`, and `tables/TABLE_H21_tds_validation.csv`. The candidate TOML and SHA sidecar remain unchanged.

## Blinded comparison

The preblind ExpH predictor retained 7.380 MW and was not used as a seed from prior experiments. ExpG's frozen robust candidate retained 3821.771 MW and converted 1580.990 MW; ExpH's analytical robust candidate converts 5391.561 MW. The ExpH requested \(\beta=10^{-4}\) is larger than ExpG's design \(\beta=1.6991\times10^{-6}\). ExpE is shown only as the pasted nominal reference (1.189211882 MW retained): no ExpE final candidate artifact was present, so it is not an independently re-read candidate and was not used as a seed.

## Overall interpretation

ExpH supports the corrected local theory: after exact gauge quotienting, the physical zero is simple and its retained-SG authority is gain-dependent; the graph/self-energy and Woodbury structure are quantitatively useful. It does **not** support the stronger claim that the complete robust co-design is solved or independently validated. The exact local predictor needs large correction, frequency excursion remains active, gains do not materially retune, no KKT/SOSC/global certificate exists, and PowerDynamics does not reproduce the requested spectral margin. The scientifically justified status is therefore **PARTIALLY_SUPPORTED**, with an analytically feasible but independently margin-failing candidate.
