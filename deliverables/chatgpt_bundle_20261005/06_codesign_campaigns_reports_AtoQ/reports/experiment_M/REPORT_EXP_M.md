# Experiment M — PD ↔ analytical model identity audit

## Verdict

**EXP_M_STATUS: FAIL / MODEL_NOT_RECONCILED.** The unchanged ExpA–K analytical model is identical to the initialized PowerDynamics (PD) model for the all-SG case after differential/algebraic reduction, but it is **not** identical at the initialized mixed SG/GFL cases. No candidate, controller gain, or earlier artifact was changed. No optimization was run.

The primary mismatch is the operating-point allocation of P/Q between colocated SG and GFL devices. The analytical model freezes the original SG and full-power GFL internal equilibria and scales their terminal currents. PD's componentwise initialization solves a different mixed equilibrium: the SG rating scales, but its operating P/Q generally does **not** scale with the retained fraction. At ExpK bus 39, the PD SG produces **78.405 MW**, while the analytical candidate assigns it **2.710 MW**. At bus 38, the PD SG produces **−0.230 MW**, while the analytical candidate assigns **+0.598 MW**. Direct SG, GFL, ZIP, and net bus P/Q sum to within `8.8e-11 MW` and `3.9e-11 MVAr` in the three audit cases ([direct-device ledger](tables/TABLE_M04_direct_device_PQ_balance.csv)).

A second, distinct issue is that pure GFL replacements retain `busbar₊Vbase = 1.0 kV` from the compiled template where `bus.csv` specifies `16.5 kV`. The corresponding line-end bases inherit this value. System Sbase and frequency are consistently **100 MVA and 60 Hz**. The full [base ledger](tables/TABLE_M01_component_bases.csv) contains 347 rows, 27 failures including duplicated bound-device rows. In a direct ExpK ablation, setting the affected busbar Vbase values from `bus.csv` left the initialized residual, spectral abscissa, and every finite pole **unchanged to the reported precision** ([ablation](tables/TABLE_M13_voltage_base_ablation.csv)). This Vbase defect therefore does not explain the observed pole shift in the current installed model.

## Frozen source and cases

The audit uses Julia 1.11.9, PowerDynamics 5.0.0 and NetworkDynamics 1.3.0. Version, package tree hashes, repository commit, installed IEEE-39 source hash, Project/Manifest hashes, and the six input CSV hashes were frozen before model construction in [SOFTWARE_PROVENANCE.toml](SOFTWARE_PROVENANCE.toml). IEEE-39 is the example bundled with the installed PowerDynamics tree, not a separately installed IEEE39 package. The installed example sets `Sbase=100` and `fbase=60` before compiling components and restores global defaults afterward. The actual input hashes match the frozen ExpG hashes ([input ledger](tables/TABLE_M03_input_hashes.csv)).

The mandatory cases are the official all-SG system, the frozen ExpG candidate, and the frozen ExpK nominal candidate ([definitions](tables/TABLE_M05_case_definitions.csv)). The installed buses use Sauer–Pai machines, AVRTypeI and TGOV1 for controlled machines, a bus-39 uncontrolled machine, ZIPLoad, and PiLine_fault ([inventory](tables/TABLE_M02_model_inventory.csv)). Each case's initialized parameters, complete compiled component matrices, state map, descriptor matrices, and poles are exported under `matrices/<case>/`.

## Gates M1–M3

### Bases

The direct initialized parameter vectors give every bus and line `Sbase=100 MVA`, `ωbase=376.99111843077515 rad/s`, and `ωframe=1`. The stock SimpleGFLDC template and mixed SG/GFL bus inspected at bus 38 also store 60 Hz. Thus the proposed 50/60 Hz construction mismatch is **not present**. Pure replacement buses have the Vbase defect described above. The installed LineEnd source binds its Vbase to the incident bus, so affected line ends remain internally consistent with the erroneous bus base and the library's consistency check does not detect the discrepancy with `bus.csv`.

### Buses 31 and 39

The ZIP outputs were read from the initialized PD state, including their initialized Vset. They are not the nominal CSV setpoints:

| Bus | SG/generator P, MW | ZIP P, MW | Net P, MW | ZIP Vset, pu |
|---:|---:|---:|---:|---:|
| 31 | 511.718874110128 | −0.107802207330 | 511.611071902798 | 9.071760316660 |
| 39 | 271.042215868631 | −375.042215868631 | −104.000000000000 | 1.767183595366 |

The Q values and machine-precision balances are in [TABLE_M04](tables/TABLE_M04_component_PQ_sharing.csv). These initialized ZIP values agree with the existing analytical KCL reconstruction. The load at either bus is not the primary explanation for ExpK's pole shift.

### Fractional architecture

The ten requested fixed-ρ cases at buses 36 and 38 were initialized and examined at `ρ = 0, 0.01, 0.5, 0.99, 1` ([TABLE_M14](tables/TABLE_M14_fractional_rho_semantics.csv), [figure](FIG_M02_internal_poles_vs_rho.png)). At interior ρ, the SG's `Sn` is `(1−ρ)Sn₀`, its `H` remains fixed, and hence `H·Sn` scales with retained rating. The GFL has no independent device-rating parameter in this implementation; its externally injected filter current is multiplied by `PortScale=ρ`. The SG and GFL state blocks are absent at their respective zero-share endpoints.

The actual dispatch differs sharply from `(1−ρ)P₀`. At bus 38, `ρ=0.5` leaves 832.585 MW on the SG instead of the 415 MW assumed by the analytical model. At `ρ=0.99`, the SG produces −84.669 MW. These are initialized PD equilibria, with the latter violating the machine's declared nonnegative torque bounds. They demonstrate that rating and terminal-current scaling alone do not enforce the intended P/Q split.

## Linear model and spectral comparison

All 255 bus/line component linearizations returned M/A/B/C/D successfully. The complete initialized parameter dump contains 4,795 rows ([parameters](tables/TABLE_M06_component_parameter_dump.csv), [dimensions](tables/TABLE_M07_component_dimensions.csv)). State inventories match in all three cases: PD and analytical total dimensions are 192, 195, and 186; their differential dimensions are 114, 117, and 108 respectively ([inventory](tables/TABLE_M09_state_inventory_comparison.csv), [physical-state alignment](tables/TABLE_M10_state_alignment.csv)).

| Case | Reduced A relative error | Maximum matched finite-pole error, s⁻¹ | PD rightmost after one gauge, s⁻¹ | Analytical rightmost, s⁻¹ |
|---|---:|---:|---:|---:|
| all-SG | 2.919e−15 | 2.302e−12 | −0.098065409336 | −0.098065409336 |
| ExpG | 2.484e−3 | 21.029 | about −9.28e−11 | −0.137932687804 |
| ExpK | 2.511e−4 | 5.279 | +0.005106001794 | −0.050001000474 |

The maximum pole errors are optimal one-to-one assignments across the complete finite spectra; they should not be interpreted as an individual branch sensitivity ([pole ledgers](tables/TABLE_M16_pole_matching_summary.csv)). The old `stability_audit` excludes **every** `|λ|≤1e−8` as gauge. In ExpG, PD has two such poles, while the all-SG case has one. Treating both ExpG poles as gauge creates an apparent rightmost-pole agreement near −0.138; only one gauge was removed for the table above. The classification of the second near-zero pole remains unresolved and is not silently discarded.

The all-SG *reduced* A and the sampled bus-port transfers match to numerical precision. The raw all-SG descriptor A matrices differ because the analytical reconstruction and compiled PD use different algebraic representations; no arbitrary dense similarity transform was used to erase this difference. See [matrix identity](tables/TABLE_M11_matrix_identity.csv), [block ledger](tables/TABLE_M12_block_error_ledger.csv), and [heatmap](FIG_M01_A_error_heatmap.png).

The official **closed-network** PD linearization is autonomous and supplies `M,A` without a defined external `B,C,D` interface. The old analytical `AN_B_port/AN_C_port/AN_D_port` files describe an **open bus-port** interface. TABLE_M11 records their direct subtraction as `NOT_COMPARABLE` instead of reporting a meaningless norm. Comparable port transfers are tested separately below; an independently defined load-input `B` is tested in the trajectory audit.

At 51 logarithmically spaced frequencies from `1e−4` to `1e3 rad/s`, all-SG bus ports 31, 36, 38, and 39 agree to at worst `5.7e−13` relative error. Pure GFL bus 31 in ExpK agrees to `4.7e−12`. The mixed bus 36 in ExpG differs by up to `0.698`, and mixed buses 38 and 39 in ExpK by up to `0.00954` and `0.793`, respectively ([port ledger](tables/TABLE_M15_port_transfer_identity.csv), [frequency plot](FIG_M03_GFL_port_error_vs_frequency.png)). This locates the mismatch in the mixed operating-point-dependent ports.

The ExpK PD pole at `+0.005106001794 s⁻¹` is present in the direct descriptor and in the independent PD-port closure. Its participation is approximately 49% GFL filter, 27% GFL PLL, and 24% SG machine, with the bus-39 SG angle the largest individual state ([TABLE_M17](tables/TABLE_M17_critical_pole_participation.csv)). Global pole matching pairs it with an analytical pole at `−0.10413 s⁻¹`; the nearest real analytical pole is `−0.05000 s⁻¹`. A unique branch correspondence has not been certified. The plain-language cause is **different mixed SG/GFL initialized dispatch and internal states**, most pronounced at bus 39.

## Corrected local analytical reconstruction

The installed `open_loop_linearization` gives exact bus-port `Zbus` and line/network `Ynw`. Its installed tests use **positive** feedback for this sign convention. Applying the positive-feedback Schur formula produces a separate `ANALYTICAL_PD_EXACT` descriptor with M/A/B/C/D and an explicit bus-group state permutation in each case directory. It reproduces the direct PD M exactly; the relative A errors are `3.34e−17`, `2.37e−20`, and `7.69e−18` for all-SG, ExpG, and ExpK. Maximum complete finite-pole differences are `7.96e−13`, `7.66e−11`, and `1.65e−9 s⁻¹` ([corrected identity](tables/TABLE_M21_corrected_closure_identity.csv), [open-loop ledger](tables/TABLE_M08_PD_open_loop_closure_identity.csv)).

This corrected reconstruction is **local to each initialized PD operating point**. It is not yet a parameterized model of how dispatch and states should change with `(ρ,Kp,Ki)`. It therefore cannot validate the old ExpK optimum or authorize another optimization.

At ExpK, the corrected closure and direct PD linearization agree in a tiny bus-16 load-step trajectory to relative `1.38e−14` and a `1e−7 pu` SG-39 speed perturbation to `2.01e−12` ([linear trajectory ledger](tables/TABLE_M19_linear_trajectory_identity.csv)). The load-input B was taken from the direct PD parameter Jacobian and shared by the two compared state models. Independent full-network current/voltage input channels were not exported. A second ForwardDiff pass agrees with the official linearization's A to `4.66e−10` maximum absolute element error ([AD check](matrices/ExpK_nominal/PD_direct_AD_check.csv)).

The independent nonlinear bus-16 load pulses at fractions `1e−3`, `1e−4`, and `1e−5` agree with the direct PD linear model to maximum relative frequency errors `1.11e−4`, `1.11e−5`, and `8.35e−6`, respectively ([scaling ledger](tables/TABLE_M20_nonlinear_linear_scaling.csv), [traces](FIG_M04_error_scaling.png)). The upper two absolute errors differ by a factor of approximately 100 when the pulse differs by 10, consistent with quadratic small-signal behavior there. A three-point log fit gives exponent `1.56`; the smallest pulse reaches a numerical error floor, so this is not a clean asymptotic exponent estimate ([fit](tables/TABLE_M20_scaling_fit.csv)). These traces validate the direct PD linearization locally, while the original mixed analytical model remains nonidentical.

## Acceptance and next action

**MODEL_IDENTITY_CERTIFIED: NO. MODEL_SAFE_FOR_Z_OPTIMIZATION: NO.** The old analytical model fails mixed matrix, port, and complete-spectrum identity. A numerically exact fixed-point PD-port reconstruction now exists, but the model generator must first impose the intended per-device P/Q sharing and reject SG equilibria that violate declared bounds. The pure GFL bus voltage base must also be set from `bus.csv` for correct physical units, although its ablation has no spectral effect here. Those changes belong in a new reconstruction; ExpA–K artifacts and frozen candidates remain untouched. After the corrected parameterized model matches PD over support and ρ sweeps, optimization may be reconsidered.

No grey-box parameter estimation was used because all needed installed parameters are exposed. No push was performed.
