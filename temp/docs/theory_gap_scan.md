# Theory Gap Scan: Master Thesis-Level Directions

Project: spectral damping-margin analysis for IBR-dominated power networks.

Purpose: review the supplied documents and extract mathematically defensible
research gaps that could extend the current C5/C7 paper into a master's
thesis-level contribution. This memo follows `HANDOFF_Codex.md`: no novelty
claim is final until targeted literature search and benchmark validation.

---

## Documents Reviewed

Local PDFs:

- Guanghui Wen et al., *Cooperative Control of Complex Network Systems with
  Dynamic Topologies* (2021).
- Lei Ding, Qing-Long Han, Boda Ning, *Distributed Control and Optimization of
  Networked Microgrids: A Multi-Agent System Based Approach* (2022).
- Zhongkui Li, Zhisheng Duan, *Cooperative Control of Multi-Agent Systems: A
  Consensus Region Approach* (2014).
- `Master_Plan_Paper20pp.pdf`.

External literature checked for novelty risk:

- System-theoretic low-inertia performance metrics:
  https://arxiv.org/abs/1703.02646
- Heterogeneous performance metrics and inertia placement:
  https://arxiv.org/abs/1710.07195
- Optimal placement of inertia and primary control via perturbation:
  https://arxiv.org/abs/1906.06922
- Full inverter-network small-signal / singular perturbation modeling:
  https://bbjohnson.org/wp-content/uploads/2022/04/Singular_Perturbation_and_Small-signal_Stability_for_Inverter_Networks.pdf
- Load/network dynamics in GFM stability:
  https://arxiv.org/abs/2212.08147
- GFL/GFM gain, line-dynamics, and operating-condition effects:
  https://arxiv.org/abs/2311.12152
- Delay-robust distributed frequency control:
  https://arxiv.org/abs/2307.14098
- Closed-form inverter-integrated small-signal stability condition:
  https://arxiv.org/abs/2503.20276
- Stability manifolds for converter-dominated systems:
  https://arxiv.org/abs/2605.26254
- GFL/GFM switching and small-signal security regions:
  https://arxiv.org/abs/2603.19618

---

## High-Level Reading

The supplied books do not replace the current C5/C7 direction; they suggest
three mathematically mature ways to extend it:

1. Li-Duan consensus-region theory gives a clean template:

   ```text
   network stability = all nonzero graph eigenvalues lie inside a region
   determined by local agent dynamics and feedback/control parameters.
   ```

   For this project, the analogous object is not a consensus region for
   first-order agents, but a damping-margin region for the reduced QEP.

2. Wen's dynamic-topology theory gives the right tools for topology changes:
   switching graphs, common/multiple Lyapunov functions, dwell-time logic,
   directed spanning-tree conditions, and resilience under attacks.

   For this project, the relevant transfer is not generic MAS consensus. It is
   stability of a family of physical grid graphs under edge switching,
   contingency sequences, or GFL/GFM mode switching.

3. Ding-Han-Ning's microgrid book gives power-system-specific extensions:
   communication delays, event-triggered updates, secondary control,
   finite-time restoration, and DoS/resilience.

   For this project, the most useful transfer is a reduced-order
   delay-dependent damping margin, not a full secondary-control redesign.

---

## Ranked Gap Candidates

### G1. Damping Consensus Region for Reduced QEPs

Core idea:

Convert Li-Duan's consensus-region idea into a power-system damping-region
criterion. For linear MAS, one asks whether `c lambda_i(L)` lies inside a
stability/performance region. Here, define a reduced-QEP damping region:

```text
R_zeta = {alpha: zeta_min(s^2 I + s D_r(alpha) + L_r(alpha)) >= zeta_star}
```

where `alpha` may represent line scaling, IBR penetration, PLL damping, virtual
inertia, or a GFL/GFM mode parameter.

Why this is promising:

- It directly extends C5 rather than replacing it.
- It turns C5 from an estimator into a design/certification geometry.
- It connects to consensus-region theory without pretending power grids are
  simple MAS consensus systems.

Mathematical route:

1. Use C5 to select the low-damping subspace `V`.
2. Parametrize topology/control changes:

   ```text
   L(alpha) = L0 + sum_e alpha_e B_e
   D(alpha) = D0 + sum_i beta_i E_i
   M(alpha) = M0 + sum_i mu_i F_i
   ```

3. Build:

   ```text
   L_r(alpha) = V^T M(alpha)^(-1/2) L(alpha) M(alpha)^(-1/2) V
   D_r(alpha) = V^T M(alpha)^(-1/2) D(alpha) M(alpha)^(-1/2) V
   ```

4. Define feasible region by:

   ```text
   min_j -Re(s_j(alpha)) / |s_j(alpha)| >= zeta_star
   ```

Controls to run:

- Compare against full ANDES eigenvalues, not only synthetic QEPs.
- Compare against low-frequency mode selection and homogeneous damping formula.
- Hold operating point fixed when isolating topology/control effects.

What is not tested:

- Whether the reduced region has stable shape under high IBR penetration.
- Whether one fixed selected subspace `V` remains valid across a large parameter
  sweep.

Originality status:

- Not confirmed. Stability/security regions for IBRs are active in 2025-2026.
  The possible gap is the reduced spectral-QEP construction and its planning use,
  not the general idea of a stability region.

Confidence:

- High as a master's thesis direction.
- Medium as a novelty claim until IEEE Xplore/arXiv search is completed.

---

### G2. Dynamic Hidden Margin via Reduced H-infinity Damping Region

Core idea:

Li-Duan's H-infinity consensus region suggests a disciplined version of the
"hidden margin" idea. Do not define hidden margin as OPF or steady-state
reserve. Define it dynamically:

```text
gamma_r = sup_omega sigma_max(
  C_r (-omega^2 I + j omega D_r + L_r)^(-1) B_r
)
```

Then a grid can have acceptable eigenvalue damping but still have a large
directional disturbance gain for certain bus/line disturbance shapes.

Why this is promising:

- It answers the user's question: "what perturbations excite the network?"
- It uses a mature H-infinity framework instead of informal resonance language.
- It can be built on C5's reduced subspace, so it remains computationally cheap.

Mechanism to test:

- The bad perturbation must come from modal alignment/input-output geometry,
  not merely from lower damping ratio.
- The correct control is equal `zeta_min` or fixed eigenvalue spectrum with
  different `B,C` disturbance/output maps.

Controls to run:

- Equalize damping ratio before comparing gains.
- Compare `gamma_r` against full state-space frequency response.
- Compare bus injection, line fault, and load-step `B` maps separately.
- Check whether the worst singular vector corresponds to a physical disturbance.

What is not tested:

- Whether the reduced resolvent remains accurate for full IBR models with
  high-order PLL/current-control states.
- Whether ANDES can expose the right `B,C` channels without extra model work.

Originality status:

- H-infinity consensus regions are known.
- Power-system disturbance amplification metrics are known.
- The candidate gap is a C5-compatible reduced-QEP hidden dynamic margin for
  IBR damping planning. Originality not confirmed.

Confidence:

- High as an extension figure/section.
- Medium as a standalone theoretical contribution.

---

### G3. Switched-Topology Damping Margin with Dwell-Time Certification

Core idea:

Wen's dynamic-topology theory suggests replacing static validation with a family
of graphs:

```text
G_sigma(t) in {G0, G1, ..., Gm}
```

where each graph corresponds to line outage, line reinforcement, topology
reconfiguration, or IBR communication-mode change. The target is not consensus;
it is preservation of damping margin under switching.

Possible theorem shape:

If every reduced subsystem

```text
A_r^q = [[0, I], [-L_r^q, -D_r^q]]
```

has `zeta_min(A_r^q) >= zeta_star`, and either:

- there exists a common quadratic Lyapunov function `P > 0`, or
- multiple Lyapunov functions satisfy an average dwell-time condition,

then the switching family preserves a certified dynamic margin.

Why this is promising:

- It uses Wen's theory directly.
- It converts "topology changes that improve robustness" into a rigorous
  switching-stability problem.
- It stays inside dynamic stability, avoiding the refuted OPF/hidden-margin
  framing.

Controls to run:

- Compare with full switched simulations, not only frozen eigenvalues.
- Test whether each frozen graph is stable but switching can still degrade
  transients.
- Equalize damping before claiming a topology mechanism.

What is not tested:

- Feasibility of common Lyapunov functions on realistic 39-bus reduced models.
- Whether switching time scales in real grid operations justify dwell-time
  assumptions.

Originality status:

- Switched-system and MAS dynamic-topology theory is known.
- Small-signal security regions and GFL/GFM switching are active topics.
- The possible contribution is reduced-QEP dwell-time certification for IBR
  damping-margin screening. Originality not confirmed.

Confidence:

- Medium-high as a master's thesis chapter.
- Medium as paper contribution because proof/validation burden is larger.

---

### G4. Delay-Dependent Spectral Damping Margin

Core idea:

Ding-Han-Ning emphasize communication delays and Lyapunov-Krasovskii/LMI tools.
The transfer is:

```text
M theta_ddot(t) + D theta_dot(t) + L theta(t)
  + K_tau theta(t - tau) = p(t)
```

or, more realistically, delay enters secondary/PLL/GFM control channels. Reduce
the system using C5 and compute a delay margin:

```text
tau_max = max tau such that zeta_min(reduced delayed system) >= zeta_star
```

Possible reduced LMI:

```text
dot x_r(t) = A0_r x_r(t) + A1_r x_r(t - tau)
```

Use delay-dependent LMIs or frequency-domain root crossing on the reduced
system.

Why this is promising:

- It brings advanced math from networked control into the IBR paper.
- It is naturally relevant to microgrids and distributed IBR control.
- It gives a tangible planning quantity: allowable delay before damping margin
  collapses.

Controls to run:

- Compare against full delayed state-space or DDE simulations where possible.
- Separate physical damping loss from delay-induced phase lag.
- Test multiple delay locations: PLL measurement, secondary control,
  communication graph, and protection signal.

What is not tested:

- Whether the main ANDES models expose delay states easily.
- Whether delay enters the primary IBR dynamics or only secondary controls in
  the intended benchmark.

Originality status:

- Delay-dependent microgrid stability is well published.
- The candidate gap is reduced spectral damping-margin estimation for IBR
  planning, not delay stability itself.

Confidence:

- Medium. Strong mathematically, but scope may become too large for the first
  paper.

---

### G5. Resilience Margin Under Edge/Node Attacks

Core idea:

Wen and Ding both treat malicious attacks and topology loss. Transfer this to
power-grid dynamics by defining an attack/contingency margin:

```text
R_k = min_{removal set S, |S| <= k} zeta_min(G \ S)
```

or the reduced approximation:

```text
R_k^C5 = min_{S, |S| <= k} zeta_min_C5(G \ S)
```

Use C7 to avoid brute force by ranking edge removals or reinforcements per
critical mode.

Why this is promising:

- It gives "weak links" a precise meaning.
- It connects resilience/protection language to dynamic stability.
- It can generate strong figures: best/worst edge removals, margin collapse,
  and targeted reinforcement.

Controls to run:

- Compare edge removal against random removal and against `lambda_2` ranking.
- Distinguish physical line removal from communication-link DoS.
- Confirm margin degradation is not just islanding/connectivity failure.

What is not tested:

- Multi-edge interactions.
- Whether the C7 first-order map predicts finite removals accurately.

Originality status:

- Resilience and attack models are known.
- The possible gap is damping-margin-targeted spectral attack/reinforcement
  screening for IBR transition planning. Originality not confirmed.

Confidence:

- Medium-high as an application section.
- Medium as theory unless a submodularity/guarantee result is added.

---

### G6. Matrix-Weighted / dq Laplacian Theory for IBR Networks

Core idea:

The scalar Laplacian `L` is a useful swing-equation abstraction, but full
inverter networks naturally produce matrix-weighted admittance/Laplacian
objects in dq/DQ coordinates. A master's-level theoretical extension is:

```text
L_scalar -> L_dq = (B kron I2) A_dq (B kron I2)^T
```

where edge weights are 2x2 admittance matrices, not scalars. Then define a
matrix-weighted reduced QEP or first-order reduction for IBR modes.

Why this is promising:

- It addresses a real limitation of the current scalar spectral model.
- It aligns with recent warnings that line/network dynamics matter for
  inverter-dominated stability.
- It can explain when the scalar C5 approximation is valid or fails.

Controls to run:

- Compare scalar lossless Laplacian, algebraic network model, and dynamic-line
  dq model on the same operating point.
- Hold inverter gains fixed when changing line model fidelity.
- Check whether C5 selected modes still predict the worst full mode.

What is not tested:

- A clean reduced-QEP form for full matrix-weighted inverter dynamics.
- Whether ANDES extraction gives enough dq network information for this route.

Originality status:

- Matrix-weighted inverter network modeling is known.
- The gap is a bridge from scalar spectral damping margin to matrix-weighted
  IBR spectral reduction. Originality not confirmed.

Confidence:

- High scientific value, high effort.
- Better as a thesis extension than as the first submitted paper.

---

### G7. Event-Triggered Stability Screening

Core idea:

Use event-triggered control logic from Ding, but apply it to monitoring/planning:
do not recompute full eigenspectra continuously. Trigger C5 updates only when
a low-cost proxy changes enough:

```text
||Delta L_tilde||, ||Delta D_tilde||, |Delta rho|, or predicted Delta zeta
```

Potential trigger:

```text
recompute C5 if |q_k^T Delta D_tilde q_k|/(2 sqrt(nu_k))
              + |q_k^T Delta L_tilde q_k|/(4 nu_k^(3/2)) > epsilon
```

Why this is promising:

- It converts event-triggering into an eigen-analysis resource-saving tool.
- It is practical for operations/planning dashboards.
- It complements C5's speedup story.

Controls to run:

- Compare against fixed-period full eigensolve and fixed-period C5.
- Measure missed events: cases where trigger stayed silent but true margin fell.
- Use pre-registered thresholds.

What is not tested:

- Realistic streaming data availability.
- Whether this belongs in a dynamics paper or an engineering/tooling paper.

Originality status:

- Event-triggered microgrid control is known.
- Event-triggered damping-margin screening may be less explored, but this is not
  confirmed.

Confidence:

- Medium. Useful, but less fundamental than G1/G2/G3.

---

## Recommended Thesis Architecture

The strongest master's thesis path is:

1. Core paper: C5 + C5 conversion ranking + C7 on IEEE 39-bus ANDES.
2. Theory extension: G1 reduced damping-region framework.
3. One advanced chapter:
   - choose G2 if the goal is disturbance/resonance theory;
   - choose G3 if the goal is topology/reconfiguration theory;
   - choose G4 if the goal is delay/networked-control theory.

Avoid trying to include all extensions in the first paper. The first paper should
remain tight:

```text
non-proportional damping -> damping-selected reduced QEP -> planning decisions
```

Then use this memo as the thesis expansion map.

---

## Immediate Experiments to Decide the Best Gap

Experiment 1: C5 region sanity

- Sweep one line weight, one damping parameter, one inertia parameter.
- Plot true `zeta_min` and C5 `zeta_min`.
- Check whether the reduced feasible region matches the full region.

Experiment 2: H-infinity hidden dynamic margin

- Construct reduced resolvent.
- Compare worst singular gain against full linear model.
- Repeat at equal `zeta_min` to rule out low-damping-only explanations.

Experiment 3: switching topology

- Choose 3-5 plausible topology states.
- Check frozen eigenvalues, common Lyapunov feasibility, and switched simulation.
- Compare against random switching and dwell-time variations.

Experiment 4: delay margin

- Add a delay channel in the reduced model.
- Estimate `tau_max` by reduced root crossing or LMI.
- Validate against full delayed simulation if feasible.

Experiment 5: resilience removal

- Remove one edge at a time.
- Compare C7-predicted weak links, `lambda_2` weak links, and true
  `Delta zeta_min`.
- Then test two-edge removals to expose interaction failures.

---

## Decision Table

| Candidate | Value | Risk | Fit with C5/C7 | Recommendation |
| --- | --- | --- | --- | --- |
| G1 damping region | High | Medium novelty risk | Excellent | Best thesis theory |
| G2 H-infinity hidden margin | High | Medium mechanism risk | Excellent | Best resonance/disturbance extension |
| G3 switched topology | High | Proof burden | Good | Strong if topology is the theme |
| G4 delay margin | Medium-high | Scope burden | Good | Use as future work unless time allows |
| G5 resilience removal | Medium-high | Can become application-only | Excellent | Good application chapter |
| G6 matrix-weighted dq Laplacian | Very high | High effort | Medium | Ambitious thesis extension |
| G7 event-triggered screening | Medium | Tooling paper risk | Good | Secondary engineering contribution |

---

## Bottom Line

The best new-theory direction is not to abandon C5/C7. It is to reinterpret C5
through the consensus-region lens:

```text
Classical MAS:
  local agent dynamics + graph eigenvalues -> consensus region

This thesis:
  inverter damping/control + inertial graph modes -> damping-margin region
```

Then add one of:

- H-infinity reduced hidden dynamic margin,
- switched-topology dwell-time damping margin,
- delay-dependent reduced damping margin.

None of these is confirmed original yet. The most honest claim today is:

> The supplied networked-control theory suggests several plausible master
> thesis-level extensions of the validated C5/C7 framework. The strongest is a
> reduced-QEP damping-region theory for IBR planning, but its originality and
> accuracy must be checked against recent stability-region work and against
> full ANDES IEEE 39-bus eigenanalysis.

