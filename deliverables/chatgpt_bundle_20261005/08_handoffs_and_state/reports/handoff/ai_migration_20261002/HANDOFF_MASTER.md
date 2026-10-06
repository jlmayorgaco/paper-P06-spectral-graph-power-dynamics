# IEEE-Style Research Handoff: SG-to-GFL Replacement, PLL Co-Design, and Beyond Nodal Damping

**Internal project-transfer dossier.** Prepared 2026-10-02 for research continuity across AI systems. This is not a submitted IEEE manuscript and does not assert established novelty or global optimality.

**Repository snapshot:** branch research/expQ2B-secure-optimum; HEAD d0fecb3264aeb855ab6700fc6ef66120d4b6c33c. The working tree contains extensive tracked modifications and untracked research artifacts. Do not clean, reset, checkout, or overwrite existing result folders. Read this handoff before running experiments.

## Abstract

This dossier consolidates the repository and decisions from this conversation around spectral graph power dynamics, “Beyond Nodal Damping,” and synchronous-generator-to-grid-following-inverter (SG-to-GFL) replacement with PLL co-design. It separates workstreams whose percentages and guarantees are not interchangeable: a thesis program on graph-coupled dynamics; TX3, TX4, and contextual-dynamic-weakness manuscripts; the frozen IAS26-110 poster package; analytical and PowerDynamics replacement experiments; and the new fixed-mix PLL counterfactual performed for this conversation.

The newest physical counterfactual fixes the ten-bus IEEE-39 replacement vector at 90.7854% GFL and compares co-designed PLL gains with a common nominal PLL. Across six load steps, retuning lowers every tested peak frequency deviation and changes frequency-threshold failures from 2/6 to 0/6, while raising RoCoF in all six events. Worst frequency deviation changes from 0.516675 to 0.499003 Hz; worst RoCoF changes from 0.242987 to 0.334563 Hz/s. Voltage passes the exploratory range in both arms. This is finite evidence that PLL tuning changes performance at fixed replacement, not proof that tuning increases maximum feasible replacement.

The mathematical direction is a network-aware finite compensability boundary: maximize weighted GFL replacement subject to declared DAE stability and operational constraints, with PLL gains as recourse variables and a declared perturbation/uncertainty set. DAE sensitivities and a first-order recourse envelope are implemented; an analytic modal derivative was checked against centered finite differences. A second-order graph-structured curvature bound and useful finite infeasibility certificate have not been established. The local tangent predicts only 0.001696 percentage points of uniform additional replacement and is not a physical capacity bound. The dossier records failures constraining stronger claims, including a prior analytical/PowerDynamics model mismatch, a graph-diagonalization gate failure, endpoint mode/dimension issues, ambiguous frequency metrics, incomplete robustness/global-search certificates, and unmodeled explicit PLL latency.

**Index terms:** power-system dynamics, grid-following converter, PLL, synchronous-generator replacement, IEEE-39, graph Laplacian, DAE sensitivity, robust optimization, project handoff.

## 1. Fast restart: claims that are safe today

The immediate research objective is to determine how much SG dispatch can be replaced by GFL resources when each installed GFL receives nodewise PLL gains, subject to network-wide dynamic and physical constraints. The theory remains parameterized by the admissible disturbance class, operating uncertainty, horizon, and hardware limits. The six bus/load combinations used so far are validation scenarios only; they do not define the general theory.

| Quantity | Common nominal PLL | Co-designed PLL | Interpretation |
|---|---:|---:|---|
| Fixed GFL share | 90.785446% | 90.785446% | Same ten-bus replacement vector |
| Worst tested max \(|\Delta f|\) | 0.516675 Hz | 0.499003 Hz | Lower for tuned PLL in all 6 events |
| Frequency failures at 0.5 Hz | 2/6 | 0/6 | Exploratory threshold, not a grid-code claim |
| Worst tested RoCoF | 0.242987 Hz/s | 0.334563 Hz/s | Tuning raises maximum by 37.7% |
| Voltage range | 0.9-1.1 pu passed | 0.9-1.1 pu passed | Only these six runs |
| Explicit hard current/DC limits | Not certified | Not certified | Traces recorded; limits not declared |

Safe wording: “The tuned PLL improves peak frequency deviation in six specified IEEE-39 load-step tests at a fixed replacement mix, with a RoCoF tradeoff.” Do not say “more stable overall,” “robust to any disturbance,” “maximum replacement,” “optimal,” or “globally optimal.”

A separate earlier iteration reported 90.0470% GFL with joint replacement/tuning versus 89.8604% with fixed gains, a 10.085 MW difference at the stated iteration budget. Four of five external event tests failed. This is a distinct run, candidate, and contract. Do not pool its percentages with the 90.7854% candidate or the 91.885% best-found ExpQ2B point.

## 2. Mathematical formulation and current theory boundary

Let \(x\) collect dynamic states, \(z\) algebraic network variables, \(\rho_i\in[0,1]\) the SG-to-GFL replacement fraction at site \(i\), \(\kappa_i=(\log K_{p,i},\log K_{i,i})\), \(w\) an exogenous disturbance, and \(\theta\) uncertain parameters. The model is

\[
\dot{x}=f(x,z,\rho,\kappa,w,\theta),\qquad 0=g(x,z,\rho,\kappa,w,\theta).
\]

For a regular index-1 algebraic branch, \(g_z\) is nonsingular and locally \(z=\phi(x,\eta,w,\theta)\), where \(\eta=(\rho,\kappa)\). The reduced dynamics are \(F(x,\eta,w,\theta)=f(x,\phi(x,\eta,w,\theta),\eta,w,\theta)\). The exact local Jacobian and parameter forcing are

\[
A=F_x=f_x-f_zg_z^{-1}g_x,\qquad
B_\eta=F_\eta=f_\eta-f_zg_z^{-1}g_\eta.
\]

First-order trajectory sensitivities satisfy \(\dot S_\eta=A S_\eta+B_\eta\). For a smooth active output branch \(h_j\), stack constraints as \(h(\rho,\kappa,w,\theta)\le0\). A local PLL recourse calculation is

\[
\max_{\tau,\Delta\kappa}\ \tau,\quad
h+J_\rho d\,\tau+J_\kappa\Delta\kappa\le0,\quad
\Delta\kappa\in\mathcal K.
\]

This is local linear feasibility, not nonlinear capacity or global search. In the current IEEE-39 calculation, \(d\) is the uniform direction, \(\mathcal K\) is a symmetric log-gain box with factor 1.1, and the tangent gives \(\tau=1.6964305\times10^{-5}\), or 0.001696 percentage points. The tangent is near an active modal constraint; trajectory rows use an SDIRK tangent. No uniform second-order remainder was computed.

A candidate finite exclusion certificate uses nonnegative \(\lambda\), \(q(\eta)=\lambda^\top h(\eta)\), and a verified curvature bound \(||\nabla^2q(\eta)||_2\le L_\lambda\) on an entire compact box \(\mathcal B\). Then

\[
q(\eta_0+\delta)\ge q(\eta_0)+\nabla q(\eta_0)^\top\delta-\tfrac12L_\lambda||\delta||^2.
\]

If the minimum of the right side over the full box is positive, no point in that box is feasible for the selected constraints. Taylor’s theorem and dual separation are standard tools. Potential contribution is a tight, reproducible, network-aware \(L_\lambda\) through the nonlinear DAE and relevant hybrid modes, with a measured finite benefit over nodal/dense alternatives. This certificate has not been implemented.

The graph-dynamic thesis uses a fuller operator of the form

\[
T(s)=s^2M+sD_0+L_B+\Pi(s),
\]

where \(L_B\) captures a declared coupling backbone and \(\Pi(s)\) contains dynamic ports and controller interactions. A scalar \(\lambda_2(L_B)\) describes static connectivity; it does not encode heterogeneous damping, inertia placement, operating-point residues, PLL dynamics, non-normal amplification, or finite disturbance response by itself. A graph basis is a coordinate change unless truncation is justified by an error bound. The preregistered F2 result found that the physical conductance graph does not diagonalize the detailed reactive-coupled GFL/network pencil; its off-diagonal reactive ratio was 0.2914, above the 0.05 modal-separation gate. No universal \(K_p(\lambda),K_i(\lambda)\) law is validated.

A future delay extension must distinguish pure digital latency \(T_i\), whose linearized feedback contains \(e^{-sT_i}\), from a first-order filter \(\tau_i\dot y_i=u_i-y_i\), whose transfer is \(1/(1+s\tau_i)\). They are not interchangeable. The current PowerDynamics counterfactual has no explicit PLL delay. Its 0.005 s output step is not a measured converter sampling period or physical delay.

The general capacity objective should be written only after the operating contract is fixed:

\[
C(\rho)=\frac{\sum_i P_i^{base}\rho_i}{\sum_iP_i^{base}},\qquad
\max_{\rho,\kappa} C(\rho)
\]

subject to \(H_j(\rho,\kappa)\le0\), where \(H_j\) is a supremum over declared disturbance set \(W\), uncertainty set \(\Theta\), initial-state set \(X_0\), and horizon. Constraints may cover stability, frequency, RoCoF, voltage, current, DC energy, actuator reserves, and recovery. A finite list of tests is not a supremum guarantee. Finite-horizon safety is not automatically asymptotic or regional stability. Limits absent from the model cannot be claimed as satisfied.

## 3. Repository workstreams: keep them distinct

| Workstream | Canonical location | Status and use |
|---|---|---|
| Spectral graph / Beyond Nodal Damping thesis | reports/thesis/ | Research framework and chapters. Read its dynamic-operator and limitations sections; do not reduce the contribution to “new Laplacian.” |
| TX3 connected intervention calculus | reports/papers/tx3_connected_intervention_calculus/ | Final IEEE Transactions manuscript and claim freeze. Supported C1/C2/D1/D2 and conditional D3; nominal and load-conditioned materiality claims rejected; weak-grid materiality unresolved. |
| TX4 policy-dependent incompatibility | reports/papers/tx4_policy_dependent_incompatibility/ | Hardened 12-page IEEEtran draft. Three claims only: policy-dependent minimal incompatibility, network-closure anatomy, and actionable boundary motion. Several tempting claims are explicitly excluded. |
| Contextual dynamic weakness (CDW) | reports/papers/cdw_contextual_dynamic_weakness/ | Separate IEEE-39 plus IEEE-68 research line on contextual reinforcement ranking. IEEE-68 replication changed the headline ranking claim. Do not transfer claims to PLL SG-to-GFL co-design. |
| Frozen IAS26 poster | reports/poster/ias2026/deliverables/IAS26-110/ | Existing poster is about H4 portfolio incompatibility/closure and a synthetic operating ensemble, not the new PLL compensability experiment. IAS26-120 audits it with limitations. |
| New compensability proof of concept | experiments/compensability_ieee39_20261002/ and reports/poster/ias2026/compensability_ieee39_20261002/ | New same-mix nonlinear PLL counterfactual and local tangent diagnostic. Preliminary only. |
| Analytical / exact-model replacement program | experiments/bnd_expA through experiments/bnd_expQ2B; reports/experiment_A through reports/experiment_Q2B | Long exploratory sequence. Final summaries are more authoritative than old candidate percentages. Model identity and physical validation changed during this sequence. |
| Robust nonlinear block certificate | experiments/nonlinear_codesign_20261001/ and reports/nonlinear_codesign_20261001/ | General DAE certificate prototype; its certified input radius is microscopic in physical normalization, so it does not support operational capacity claims. |
| Iterative nonlinear design | experiments/analytic_iteration_20261001/ and reports/analytic_iteration_20261001/ | Validated local-gradient iteration and PowerDynamics checks. Joint tuning improved the training objective, but 4/5 external events failed. |
| Poster concept and novelty review | reports/poster/ias2026/concept_review_20261001/ and concept_review_20261002/ | Prior concept decisions and evidence status. IAS26-110 remains a distinct poster artifact. |

The IAS26-110 package was frozen and audited; its safe claims are in its supplement and IAS26-120 adversarial audit. It reports a nominal H4 blocker, local versus collective closure, a synthetic 967-valid-ID ensemble, fixed-retuning mixed outcomes, Python/Julia same-model parity, and conditional phasor-DAE TDS evidence. It does not claim EMT, field, universal performance, or that Python/Julia implementations are independent physical models. Keep that package intact.

## 4. Experiment chronology: what worked, what failed, what changed

### 4.1 Graph and modal branch

- F0 passed: a passive conductance port graph from the Hermitian part of the Kron-reduced branch admittance was SPD and current-equivalent.
- F1 passed for isolated commuting modes: the maximally damped scalar/modal result \(d^*=2\sqrt{\nu}\) matched numerical poles. It is a reference identity, not a heterogeneous full-network theorem.
- F2 failed its registered modal-separation gate: the conductance basis diagonalized conductance coupling, not the actual reactive-coupled GFL pencil. F3/F4 graph-gain/controller branches were stopped. Do not present the isolated modal result as the IEEE-39 PLL design law.

### 4.2 Endpoint, analytical surrogate, and model correction

- G/H produced useful endpoint and local Jordan/Feshbach diagnostics, but did not establish a full robust co-design. ExpH’s independent PowerDynamics rightmost pole failed the requested \(-0.05\,s^{-1}\) margin and differed from the analytical surrogate.
- J0 refuted a proposed square-root law for the tested one-sided near-zero branch: fitted local exponents were about 0.999 to 1.001. It also found that physical dynamic dimension jumps at \(\rho=0\), and extra SG open-loop poles can dominate the critical branch. An endpoint derivative is not a general design law.
- ExpK explored 1,024 supports and returned near-100% GFL best-found points, but ExpM later established that the old analytical model was not identical to installed PowerDynamics at mixed SG/GFL points. The optimizer was explicitly declared unsafe for reuse. Those high percentages are not valid physical capacity claims.
- ExpM found that the 50/60 Hz explanation was false: the installed example used 100 MVA and 60 Hz. It found a pure-GFL voltage-base mismatch (1.0 kV versus a specified 16.5 kV) that did not explain the pole mismatch. The important discrepancy was initialized device dispatch: SG bus 39 had 78.405 MW in PowerDynamics versus 2.710 MW in analytical sharing. An initialized-port closure reproduced PowerDynamics matrices, but was not yet a parameterized optimizer.
- ExpN established a corrected PowerDynamics-exact model stage and local fixed-support KKT diagnostics; it screened 1,024 supports and traced 40 branches, but branch completeness/globality were not certified.
- ExpP reported a fixed-support local KKT candidate and exact modal/root identities, but campaign status is FAIL_TDS. A locally certified small-signal root is not sufficient if the required nonlinear time-domain check fails.
- ExpQ/Q2 found frequency-metric dependence, non-affine steady-frequency response, and that PLL gains do not create sustained primary power response. For abrupt phase steps an unfiltered continuous RoCoF can be unbounded; filtering/window choice changes the number. Use an explicitly declared frequency estimator and window.

### 4.3 Best-found candidate versus local/global proof

ExpQ2B status is FAIL_OPTIMIZATION, despite its best-found physical candidate. It reports 91.885% GFL (4,964.315 MW) and 438.446 MW retained SG, with tested PowerDynamics event limits passing and close analytic/PD poles. But the 50,001-node robust interval bound remained incomplete, KKT stationarity residual was 0.830, SOSC was not tested, and support/global completeness remained open. Do not call it a local optimum, robust optimum, or global maximum.

### 4.4 Earlier trajectory-gradient experiment

The 2026-10-01 iterative design used \(\rho=0.875\) as its common start. Reducing trajectory-gradient timestep from 0.025 s to 0.0125 s moved relative gradient error from 2.187% to 0.5525% under an unchanged 2% gate. Joint design reached 90.0470% GFL versus 89.8604% with fixed gains. Explicit step algebra matched the QP to \(1.924\times10^{-14}\); PowerDynamics reproduced three in-sample designs. But 4/5 external events failed frequency or governor actuator margin. This shows the gradient pipeline can produce an improved local candidate, and that event generalization and actuator modeling need work.

### 4.5 New six-event PowerDynamics counterfactual

The candidate vector was fixed at 90.78544633% GFL with 497.84032003 MW SG retained. Common nominal PLL is \(K_p=2\pi\cdot5\), \(K_i=(2\pi\cdot5)^2/4\) at every GFL. The co-designed arm uses candidate-specific gains. Both arms used the same network, loads, operating point, replacement locations, and event set; each event was initialized independently in PowerDynamics and integrated for 61 s including 1 s prefault, with 0.005 s output spacing and tolerance \(10^{-9}\). The six load steps were \((8,\pm100)\), \((16,\pm100)\), and \((29,\pm100)\) MW. Theory is not restricted to these events; they are only the finite test set.

| Event | Peak \(|\Delta f|\), nominal -> tuned (Hz) | Peak RoCoF, nominal -> tuned (Hz/s) |
|---|---:|---:|
| bus 8, +100 MW | 0.447211 -> 0.429913 | 0.154497 -> 0.192318 |
| bus 8, -100 MW | 0.502778 -> 0.482508 | 0.131987 -> 0.152231 |
| bus 16, +100 MW | 0.473503 -> 0.456160 | 0.115897 -> 0.123944 |
| bus 16, -100 MW | 0.516675 -> 0.498723 | 0.118298 -> 0.122017 |
| bus 29, +100 MW | 0.464380 -> 0.450945 | 0.138748 -> 0.175471 |
| bus 29, -100 MW | 0.499597 -> 0.499003 | 0.242987 -> 0.334563 |

Tuning improves frequency peak in each row and worsens RoCoF in each row. Maximum RoCoF stays under the exploratory 0.5 Hz/s threshold in these events; this does not remove the tradeoff. All voltage samples passed the exploratory 0.9-1.1 pu range. Current ratio and DC-voltage traces were recorded in the broader physical candidate summary, but no hard device limits were certified.

![Fixed-mix PLL retuning tradeoff from six IEEE-39 events](reports/poster/ias2026/compensability_ieee39_20261002/retuning_tradeoff.png)

### 4.6 Local envelope and derivative audit

The local recourse QP uses 36 normalized rows: six modal rows plus five event metrics for each of six events. The symmetric log-gain box is factor 1.1. The uniform tangent allows \(\Delta\rho=1.6964305\times10^{-5}\), equivalent to 0.001696 percentage points. The analytic modal sensitivity and a fresh centered finite difference differ by \(2.605\times10^{-6}\) relative error; the stored updated modal gradient matches the analytic derivative. Trajectory rows remain a finite-step tangent without a second-order remainder bound. Label this output LOCAL_LINEAR_TANGENT_DIAGNOSTIC_ONLY.

## 5. Error and correction register

1. **Analytical model mistaken for installed physical model.** ExpM found dispatch/base inconsistencies and pole mismatch; the old optimizer was stopped. Start from the ExpN frozen model SHA and device/PQ contract when using that lineage.
2. **Near-100% candidates over-interpreted.** ExpK ratios predate model-identity failure; ExpP local KKT is not rescued by TDS; Q2B lacks KKT and complete interval robustness. None is a certified global physical maximum.
3. **A candidate is not a proof of the maximum.** A failed local optimizer does not prove infeasibility. A feasible local point is a lower bound on achievable capacity under its exact contract, not an upper bound. Global claims require covering supports, branches, and continuous parameter boxes or a valid upper bound.
4. **A first-order tangent is not a finite capacity bound.** The 0.001696 percentage-point envelope is not a physical limit without a verified remainder.
5. **A Laplacian basis does not make heterogeneous dynamics modal.** F2 measured 29.14% reactive off-diagonal coupling. A graph basis helps organize interactions; truncation and diagonalization claims need residual/error bounds.
6. **The all-GFL endpoint changes model dimension and has gauge modes.** Treat quotient, structural zero, open-loop poles and support changes explicitly; do not apply a one-sided root law across the dimension jump.
7. **Frequency and RoCoF depend on measurement definition.** Preserve exact signal, filter, sampling rate, window and buses. PLL frequency, SG rotor frequency, bus-voltage phase frequency, and filtered RoCoF are not interchangeable.
8. **Graph cycles, Fourier transforms, Taylor expansion, Schur complements, Farkas separation, SQP, and port-Hamiltonian notation are not novelty by themselves.** The contribution has to be the result they enable, with comparison and evidence.
9. **Delay and first-order lag are different models.** Do not fit an unexplained \(\tau\) and call it a physical digital delay. The current experiment has no explicit delay.
10. **Conversation/report corrections:** one intermediate progress note misread the bus-8, +100 MW comparison; the final event table is authoritative and shows lower frequency peak but higher RoCoF for tuning in all six. The initial Julia include had a repeated OUT name; the reproducible script now uses COMP_OUT. The summary writer initially emitted literal backslash-n separators; the Spanish report was rerun with real line breaks.

## 6. Poster line and conference state

The frozen IAS26-110 poster is an audited deliverable for its own policy-dependent incompatibility story. Its audit status is PASS_WITH_EXPLICIT_LIMITATIONS; housekeeping remained incomplete. It is one page at 36 x 48 in and includes H4 nominal portfolio evidence, closure anatomy, a synthetic operating ensemble, mixed fixed-retuning outcomes, same-model Python/Julia parity, and conditional phasor-DAE TDS. It does not show the new 90.785% PLL comparison.

The compensability experiment can become a distinct poster only after it answers more than “tuning changes the transient”: does retuning measurably increase maximum feasible SG-to-GFL dispatch under a declared multi-metric contract, and can a graph/DAE method predict or certify that capacity boundary? A plausible future title is “A Network-Aware PLL-Compensability Boundary for SG-to-GFL Replacement,” but the title is a hypothesis, not a novelty claim.

The current IAS 2026 Graduate Student Poster Competition page listed May 15, 2026 as application deadline and May 30 for decisions; the Annual Meeting is October 4-8, 2026. As of this dossier date, a new application is past the published deadline. If the author already applied/was selected, these finite results may serve as preliminary material with limits. Otherwise preserve them for a later venue. Official pages: https://ias-am.ieee.org/2026/student-poster-competition-program/ and https://ias-am.ieee.org/2026/call-for-papers/.

## 7. Priority plan for the next AI/researcher

1. **Freeze the exact contract.** Pick one model lineage, objective (dispatch-weighted GFL share), initialized SG/GFL P/Q split, frequency estimator, disturbances \(W\), uncertainty \(\Theta\), horizon, actuator/reserve assumptions, and hard current/DC limits. Keep a finite scenario screen separate from a universal/robust claim.
2. **Reproduce before extending.** Read protocol and existing CSV/TOML outputs first. Re-run the small summary/gradient audit only if needed. Do not repeat all F or TX campaigns as a starting point.
3. **Trace a real capacity frontier.** At a sequence of \(\rho\) levels, re-equilibrate and optimize local PLL gains with identical budgets/starts for nominal/fixed PLL, uniform PLL, nodewise PLL, and graph-aware method. Enforce frequency, RoCoF, voltage, mode margin, current, DC energy and actuator constraints. Include multistart results and KKT residuals; do not infer infeasibility from optimizer failure.
4. **Use independent holdouts.** Tune on one declared set and test unseen buses, directions, operating points, and uncertainty draws. Expand beyond load steps only after event models are physically defined.
5. **Prove one useful certificate on a small network.** Bound second-order trajectory remainder for an active smooth constraint branch; handle switching peak branches, algebraic regularity and limiter regions. Compare graph-walk/cycle bounds against nodal and dense bounds. If vacuous, revise or reject the hypothesis.
6. **Test delay separately.** Use implementation latency if known; otherwise define an uncertainty range and label it a screening assumption. Compare pure delay/DDE or validated Padé with a separate first-order measurement filter. Ask whether \(\rho^\star(\bar T)\) and nodewise delay sensitivities depend predictably on \(\lambda_2(L_B)\), a physically defined damping graph, or full modal participation.
7. **Cross-model validation.** Reconcile initialization and device contracts in PowerDynamics and ANDES before treating agreement as independent validation. Current Python/Julia GFL11 parity in the poster is same-model implementation parity.
8. **Review literature before novelty language.** Existing work covers multi-converter PLL retuning, grid structure/Laplacian stability on IEEE-39, PLL delay effects, and controller feasibility regions. Prospective novelty must be a demonstrated network-aware nonlinear capacity boundary/certificate, not one ingredient alone.
9. **Keep an artifact ledger.** Every quantitative sentence should point to a run, input hash, output table and gate. Preserve failed runs and preregistered amendments; never overwrite frozen evidence.

## 8. Reproduction and artifact map

### Current conversation experiment

- Physical nonlinear comparison: julia --startup-file=no --project=. experiments/compensability_ieee39_20261002/pd_compare.jl
- Local tangent envelope and modal finite-difference audit: julia --startup-file=no --project=. experiments/compensability_ieee39_20261002/linear_envelope.jl
- Summary CSV, Spanish report and tradeoff plot: C:/Users/walla/anaconda3/python.exe experiments/compensability_ieee39_20261002/summarize_results.py
- Protocol: experiments/compensability_ieee39_20261002/protocol.toml
- Outcomes: reports/poster/ias2026/compensability_ieee39_20261002/events.csv, comparison.csv, linear_envelope.toml, RESULTADOS_ES.txt, and retuning_tradeoff.png/.svg/.pdf.
- Candidate provenance: reports/nonlinear_codesign_20261001/physical_candidate_summary.json, candidate_final_physical.toml, and code under experiments/nonlinear_codesign_20261001/.
- Julia project: root Project.toml and Manifest.toml; validated ExpN environment recorded Julia 1.11.9, PowerDynamics 5.0.0, NetworkDynamics 1.3.0.

### Main prior handoff documents

- Current idea/evidence status: reports/poster/ias2026/concept_review_20261002/HIPOTESIS_Y_PRUEBAS_ES.txt and EVIDENCE_STATUS.json.
- Iterative design: reports/analytic_iteration_20261001/refined/RESULTADOS_ES.txt.
- Nonlinear certificate formulation: reports/nonlinear_codesign_20261001/coupled_dq_contract/CONTROL_MODAL_COMPARACION_Y_CODESIGN_ES.txt.
- ExpQ2B final summary and full report: reports/experiment_Q2B/FINAL_SUMMARY_EXP_Q2B.md, REPORT_EXP_Q2B.md.
- Model discrepancy: reports/experiment_M/FINAL_SUMMARY_EXP_M.md.
- Corrected model and local search: reports/experiment_N/FINAL_SUMMARY_EXP_N.md; follow-up reports/experiment_P/FINAL_SUMMARY_EXP_P.md.
- Graph modal gate: reports/experiment_F_program/FINAL_SUMMARY.md.
- Endpoint correction: reports/experiment_J0/J0_DECISION.md.
- IAS poster audit: reports/poster/ias2026/deliverables/IAS26-120/IAS26-120_ADVERSARIAL_AUDIT.md.
- TX3/TX4/CDW manuscript READMEs: reports/papers/tx3_connected_intervention_calculus/README.md, reports/papers/tx4_policy_dependent_incompatibility/README.md, and reports/papers/cdw_contextual_dynamic_weakness/README.md.
- Root structure and environment: README.md, Project.toml, Manifest.toml.

### External references verified for this dossier

[R1] L. Huang, H. Xin, W. Dong, and F. Dörfler, “Impacts of Grid Structure on PLL-Synchronization Stability of Converter-Integrated Power Systems,” IFAC-PapersOnLine, vol. 55, no. 13, pp. 264-269, 2022, doi: 10.1016/j.ifacol.2022.07.270. It uses a reduced-network eigenvalue and validates structural PLL stability conclusions on IEEE-39.

[R2] L. Guo, C. Zhao, and S. H. Low, “Graph Laplacian Spectrum and Primary Frequency Regulation,” arXiv:1803.03905, 2018. It connects graph Laplacian modes with inertia/damping and frequency response under stated assumptions.

[R3] A. Bačić et al., “A graph theoretic view on small signal stability of inverter-based power grids,” arXiv:2607.08260, 2026. Preprint; its abstract reports graph-cycle contributions typically small in three IEEE test cases under its assumptions.

[R4] “Analysis of Phase-Locked Loop Filter Delay on Transient Stability of Grid-Following Converters,” Electronics, vol. 13, no. 5, 986, 2024, doi: 10.3390/electronics13050986. It studies PLL-filter delay effects on GFL transient stability.

[R5] “Small-Signal Modeling and Stability Analysis of a Grid-Following Inverter with Inertia Emulation,” Energies, vol. 16, no. 16, 5894, 2023, doi: 10.3390/en16165894. It includes digital computation/modulation delay in an inverter model; its delay assumptions are not universal hardware values.

[R6] X. Zhou et al., “Comparative Analysis of the Power Output Capabilities of Grid-Following and Grid-Forming Inverters Considering Static, Dynamic, and Thermal Limitations,” IEEE Transactions on Power Systems, vol. 39, no. 2, pp. 2693-2705, 2024, doi: 10.1109/TPWRS.2023.3279373.

[R7] “Synchronization Stability Analysis and Parameter Design of Grid-Following Inverters Considering the Interactions of Current Control and Phase-Locked Loop,” IEEE Journal of Emerging and Selected Topics in Power Electronics, vol. 12, no. 5, pp. 5013-5027, 2024, doi: 10.1109/JESTPE.2024.3444281.

[R8] “Design of a Bandwidth Limiting PLL for Grid-Tied Inverters With Guaranteed Stability,” IEEE Transactions on Sustainable Energy, vol. 16, no. 2, pp. 774-784, 2025, doi: 10.1109/TSTE.2024.3483988.

The bibliography is a starting set, not a systematic review. Check full papers and expand the search before asserting priority or absence of prior art.
