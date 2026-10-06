# Experiment D — Single-bus SG-to-GFL replacement optimum

## Result

**EXP_D_STATUS: PASS for the preregistered single-bus objective and PLL domain.** The primary bus-33 optimum is

\[
\rho^\star=1,\qquad K_p^\star=31.4159265359\;\mathrm{rad/s},\qquad
K_i^\star=246.7401100272\;\mathrm{rad/s^2},\qquad
P_{\rm GFL,max}=632\;\mathrm{MW}.
\]

The full-replacement rightmost pole is \(-0.09884700665\pm j0.43026220071\;\mathrm{s^{-1}}\), giving a spectral margin of 0.09884700665 s⁻¹ against the required 0.05 s⁻¹. The detailed PowerDynamics pole agrees exactly at the stored precision. The same saturated optimum is feasible at buses 30, 35, and 37; see `TABLE_D18_cross_bus_optima.csv`.

## Certified design domain

The frozen preregistration specifies

\[
0\le\rho\le1,\qquad 0.9\le\beta\le1.1,\qquad
K_p=\beta K_{p0},\quad K_i=\beta^2K_{i0},
\]

with \(K_{p0}=5(2\pi)\), \(K_{i0}=K_{p0}^2/4\), and required decay rate \(\sigma_{\rm req}=0.05\;\mathrm{s^{-1}}\). The proportional and integral gains therefore move together along the preregistered physical PLL-bandwidth coordinate. Nominal gains are feasible and are the unique zero-effort controller tie-break.

## Global optimality certificate

For each bus, \(P_i^0>0\) and \(0\le\rho\le1\), so every admissible point satisfies

\[
P_i^0\rho\le P_i^0.
\]

The frozen candidate has the witness \((\rho,\beta)=(1,1)\). The full-state ExpC matrix has 111 finite modes at each full-GFL endpoint, including one gauge mode; every one of the 110 non-gauge modes satisfies \(\Re\lambda\le-0.05\). A fresh PowerDynamics equilibrium and full linearization independently verified all four witnesses. Thus the objective upper bound is attained and the global optimality gap is zero. No admissible combination can replace more than 100% of the dispatch. The lexicographic controller-effort tie-break selects \(\beta=1\), hence the reported \(K_p^\star,K_i^\star\). The objective-bound argument is exact; spectral feasibility is numerically verified from the full double-precision eigenvalue sets and independently cross-checked against PowerDynamics, with the engineering-margin decision applied using a 1e-10 s⁻¹ tolerance. No interval-arithmetic eigenvalue enclosure is claimed.

This is an exact saturated-bound certificate. We did not construct a system-level polynomial \(F(s,\rho,K_p,K_i)\), compute resultants, or enumerate nonbinding pole-boundary branches: none can improve the objective past the attained physical upper bound \(\rho=1\). Consequently, this result certifies the global maximum replacement and selected tie-break point, but does not claim a complete algebraic map of every pole-boundary branch.

## Frozen candidate and independent endpoint validation

The analytical candidate was frozen before the validation adapter imported PowerDynamics. Frozen candidate SHA-256:

```text
812a792bea2e69d3878e1347c2a00c8b372b1b55b3383ad730bd8403dfd27f01
```

PowerDynamics ran after the hash gate. All four buses pass equilibrium qualification, the 0.05 s⁻¹ all-mode margin, and rightmost-pole matching; the absolute pole and alpha errors are zero in `TABLE_D11_PD_validation.csv`. The analytic candidate computation itself uses the frozen ExpC full-state matrices and makes zero PowerDynamics calls.

## Post-freeze PowerDynamics grid

The independent validation grid contains 11 replacement shares \(\rho=0,0.1,\ldots,1\) and three preregistered PLL scales \(\beta=0.9,1,1.1\) at each of buses 30, 33, 35, and 37: 132/132 cells evaluated. Every eigenvalue for every cell is retained in `TABLE_D19_PD_grid_full_spectra.csv`; the row-level grid also stores the rightmost non-gauge pole. An audit confirms that full-spectrum mode counts and all-mode pass/fail classifications reproduce every grid row. Overall, 123 cells meet the required margin. Buses 30, 33, and 35 pass all 33 points. Bus 37 passes 24/33; all nine violations occur at \(\rho=0.7,0.8,0.9\), one for each beta, with equilibria converged. At \(\rho=1\), all three gain points pass at all four buses.

The bus-37 result is nonmonotone: the interior mixed realization has a right-half-plane or weakly damped mode at the three high-share grid columns, while the zero-share SG is structurally removed at \(\rho=1\) and the endpoint is stable. The endpoint remains the global maximum for the stated static replacement objective; the interior stability holes matter for staged replacement and must be considered in any path-dependent implementation. The figure and row-level data are `figures/FIG_D01_PD_validation_grid.png` and `tables/TABLE_D17_PD_rho_PLL_grid.csv`.

## Supporting results and limitations

- The local nine-state PLL model has exact affine dependence on its implemented gains; the Kp and Ki updates have rank one each and joint rank two. The Woodbury local-resolvent factorization and same-equilibrium local current-port update remain documented in the structural-audit artifacts.
- This ExpD optimum uses exact full-state endpoint spectra, not a graph-mode truncation or a low-rank closure. It establishes replacement capacity under the stated frozen small-signal model and controller domain.
- The grid follows the preregistered one-dimensional PLL coordinate. If the intended problem is instead an independent rectangular Kp/Ki domain, its numerical bounds must be frozen separately; the saturated objective proof is unchanged whenever the nominal gain point remains admissible, but the validation surface must be recomputed on that domain.
- The current-sharing model has no current limiter. No current-limit, nonlinear transient, EMT, uncertainty-robustness, or multi-bus co-design claim is made here.

## Decision

ExpD delivers the requested numerical optimum for bus 33 and the three cross-bus repetitions, with a zero-gap global certificate and independent endpoint PowerDynamics validation. The post-freeze grid also exposes the bus-37 interior stability hole. ExpE can begin from these endpoint results, carrying that bus-37 behavior as a constraint on any staged or partial replacement policy.
