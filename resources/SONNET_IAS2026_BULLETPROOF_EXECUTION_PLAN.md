# SONNET EXECUTION PLAN — IAS 2026 VANCOUVER “BULLETPROOF” CLOSURE CAMPAIGN

**Project:** Minimal Replacement Blockers / Robust Transition Compatibility in SG→GFL grid transitions  
**Target:** IEEE IAS Annual Meeting 2026 Vancouver — poster competition + defensible companion paper  
**Primary benchmark:** IEEE 39-bus  
**Independent implementation:** Julia / PowerDynamics.jl  
**Secondary scale benchmark:** IEEE 118-bus or another documented larger dynamic benchmark  
**Languages:** Python + Julia  
**Agent:** Claude Sonnet / Claude Code  
**Mode:** execution + audit + falsification, not brainstorming  
**Final deliverable:** one self-contained ZIP with code, environments, raw data, derived data, plots, reports, paper PDF, claim ledger, negative results, and reproducibility manifest.

---

# 0. READ THIS FIRST — NON-NEGOTIABLE OPERATING RULES

You are executing a scientific closure campaign. Do **not** optimize for producing positive results. Optimize for producing correct, reproducible, reviewer-defensible results.

## 0.1 Never modify frozen historical evidence

1. Discover the existing repository root.
2. Record:
   - current branch;
   - current HEAD;
   - all tags/branches containing `TX4`, `IAS2026`, `FINAL`, `CDW`;
   - dirty/untracked files.
3. Do **not** reset, rebase, delete, overwrite, or “clean up” historical result directories.
4. Do **not** push.
5. Create a new branch from the canonical freeze only after identifying it.

Expected historical references, to be verified rather than blindly trusted:

```text
TX4 parent freeze:
69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6

Known audit branch:
research/tx4-final-modal-scope-audit

Known CDW/topology work:
research/contextual-dynamic-weakness
```

If the expected commit does not exist locally, **do not invent it**. Record what exists and use the nearest documented canonical source only after proving provenance.

Recommended new branch name:

```text
research/ias2026-bulletproof-closure-v1
```

## 0.2 No post-hoc retuning

Once any blind/holdout prediction is frozen:

- commit the preregistration;
- save SHA-256 checksums;
- do not change thresholds, selected policies, candidate actions, or metrics after reveal;
- any deviation must be written to `docs/DEVIATIONS.md` before re-running, with reason and effect.

## 0.3 Every result must have one of these labels

Use exactly:

```text
PROVED
NUMERICALLY_VERIFIED
IEEE39_VALIDATED
POWERDYNAMICS_VALIDATED
SECOND_MODEL_VALIDATED
NONLINEAR_TDS_VALIDATED
SYNTHETIC_PILOT
CONSTRUCTED_COUNTEREXAMPLE
RETROSPECTIVE
BLIND_HOLDOUT
NOT_TESTED
REFUTED
UNRESOLVED
STOPPED_BY_GATE
```

Never call something “validated” if it only passed a synthetic toy.

## 0.4 A STOP result is a scientific result

If a gate fails:

- stop the dependent chain;
- save all evidence;
- state exactly what failed;
- do not patch around the gate merely to continue;
- independent phases that do not depend on that gate may continue.

## 0.5 Do not use “ground truth”

Use:

> full-order reference within the frozen model

unless comparing against genuinely independent physical/experimental evidence.

## 0.6 Never hide negative results

Create and continuously update:

```text
docs/NEGATIVE_RESULTS.md
```

This file is mandatory in the final ZIP.

---

# 1. SCIENTIFIC OBJECTIVE

The campaign must test the following research programme adversarially.

We study a finite family of SG→GFL replacement configurations. Let \(S\) denote the set of replaced synchronous generators and let \(\theta\) denote controller/operating policy.

The central questions are:

1. **WHO fails?**  
   Which inclusion-minimal replacement portfolios are dynamically incompatible?

2. **WHY do they fail?**  
   Is the instability caused by local device singularity, or by network-mediated collective feedback?

3. **HOW MODEL-ROBUST is the conclusion?**  
   What structured terminal-model uncertainty is required to destroy transverse stability?

4. **HOW should the transition be repaired?**  
   Can topology/control/support actions make the **entire implementation path** robustly safe, rather than only one final eigenvalue?

The campaign must attempt to falsify each layer.

---

# 2. CANONICAL THEORY TO IMPLEMENT AND AUDIT

Do not silently alter these definitions. If you find a mathematical flaw, record it explicitly and propose the corrected statement in the claim ledger.

## 2.1 Exact phasor DAE

Use the existing frozen model:

\[
\dot x=f(x,z;S,\theta),\qquad 0=g(x,z;S,\theta)
\]

with the index-one reduction

\[
A=f_x-f_zg_z^{-1}g_x.
\]

## 2.2 Transverse dynamics

The frozen no-primary-restoration model contains the rotational gauge / frequency-drift structure.

Construct and numerically verify the invariant subspace before every robust calculation.

Define the exact transverse matrix \(A_\perp\).

Required numerical residuals:

\[
\frac{\|A R_x\|}{\|A\|\|R_x\|}
\]

and, when the drift partner applies,

\[
\frac{\|Aw-\omega_B R_x\|}
{\|A\|\|w\|}.
\]

Store both.

No robust-stability calculation is allowed on a nominal operator that still contains the structural zero/Jordan chain.

## 2.3 Global and mode-scoped blockers

Full-spectrum spectral abscissa:

\[
\alpha_\perp(S,\theta)=\max_{\lambda\in\sigma(A_\perp)}\Re\lambda.
\]

Mode-scoped:

\[
\alpha_\Omega(S,\theta)=
\max_{\lambda\in\sigma(A_\perp),\,f(\lambda)\in\Omega}\Re\lambda
\]

with the explicit convention:

\[
\alpha_\Omega=-\infty
\]

if no eigenvalue exists in the band.

Use \(\Omega=[0.3,1.5]\) Hz only where this is the preregistered target family.

Define separately:

\[
\mathcal H_\perp(\theta)
=
\min_{\subseteq}\{S:\alpha_\perp(S,\theta)\ge0\}
\]

and

\[
\mathcal H_\Omega(\theta)
=
\min_{\subseteq}\{S:\alpha_\Omega(S,\theta)\ge0\}.
\]

Do not conflate them.

## 2.4 Exact terminal port operator

Use a dedicated symbol:

\[
\mathcal T_S(s)
=
g_z+g_x(sI-f_x)^{-1}f_z.
\]

Do not reuse \(\mathcal T\) for a second-order angle model.

If retained/self-energy models are used, call them \(\mathcal S_r(s)\) or \(T_\theta(s)\), never the same symbol as the exact voltage-port operator.

## 2.5 Exact replacement increment

Do not write “up to KCL sign”.

Define only:

\[
\Delta\mathcal T_i(s)
=
E_i^\top
[\mathcal T_{\{i\}}(s)-\mathcal T_0(s)]
E_i.
\]

Then:

\[
D_S=\operatorname{blkdiag}(\Delta\mathcal T_i)_{i\in S},
\qquad
K_S=E_S^\top\mathcal T_0^{-1}E_S.
\]

Split:

\[
K_S=K_{d,S}+K_{o,S}.
\]

Define

\[
L_S=I+D_SK_{d,S}.
\]

Only if \(L_S\) is nonsingular define:

\[
Q_S=L_S^{-1}D_SK_{o,S}.
\]

The local-factor nonsingularity is an explicit theorem assumption.

## 2.6 Exact closure identities

Verify:

\[
\frac{\det\mathcal T_S}{\det\mathcal T_0}
=
\det(I+D_SK_S)
\]

and:

\[
\det(I+D_SK_S)
=
\left[\prod_{i\in S}\det(I+M_{ii})\right]\det(I+Q_S).
\]

At a collective boundary, when local factors remain regular:

\[
-1\in\sigma(Q_H).
\]

## 2.7 Contextual return

For \(H=\{i\}\cup R\):

\[
R_{i|R}
=
Q_{iR}(I+Q_{RR})^{-1}Q_{Ri}.
\]

Verify:

\[
\det(I+Q_H)
=
\det(I+Q_{RR})
\det(I-R_{i|R})
\]

and therefore at the boundary:

\[
+1\in\sigma(R_{i|R}).
\]

Also verify the local critical component is nonzero:

\[
q_i\neq0
\]

for each member of a minimal blocker, under the regular proper-subset assumptions.

## 2.8 Collective–Robustness Bridge

Derive using a Woodbury identity that does **not** require \(D_S^{-1}\).

Define the portfolio terminal Green map:

\[
G_S^p
=
E_S^\top
\mathcal T_S^{-1}
E_S.
\]

Derive:

\[
G_S^p
=
K_S(I+D_SK_S)^{-1}
\]

and then:

\[
\boxed{
G_S^p
=
K_S(I+Q_S)^{-1}L_S^{-1}.
}
\]

This is the key bridge to robustness.

For a terminal uncertainty interconnection:

\[
\mathcal L_S
=
W_R\,G_S^p\,W_L
\]

obtain:

\[
\boxed{
\mathcal L_S
=
W_RK_S(I+Q_S)^{-1}L_S^{-1}W_L.
}
\]

For a fixed-frequency full-complex additive uncertainty where the reciprocal singular-value formula applies, let:

\[
A_S=W_RK_S,\qquad
C_S=L_S^{-1}W_L.
\]

Verify the bounds:

\[
\frac{\sigma_{\min}(I+Q_S)}
{\sigma_{\max}(A_S)\sigma_{\max}(C_S)}
\le r_{\rm ff}(S)
\le
\frac{\sigma_{\min}(I+Q_S)}
{\sigma_{\min}(A_S)\sigma_{\min}(C_S)}.
\]

Do **not** state “well conditioned is enough”.

The asymptotic statement requires \(A_S,C_S\) to remain:

- uniformly bounded;
- boundedly invertible.

## 2.9 Transverse robust stability

This is mandatory.

The full frozen model has a structural zero/Jordan chain. Therefore robust analysis must be done on a deflated/transverse representation.

Define:

\[
r_\Delta^\perp(S)=0
\]

for a nominally transverse-unstable portfolio.

For a nominally transverse-stable portfolio:

\[
r_\Delta^\perp(S)
=
\inf\{
\|\Delta\|:
\text{structured uncertainty creates transverse instability}
\}.
\]

Construct either:

1. an exact transverse LFT from \(A_\perp\), or
2. a symmetry-deflated exact port operator \(\overline{\mathcal T}_S\).

Do not call any \(\mu\) result valid until the structural zero has been removed and the nominal uncertainty loop is well posed.

## 2.10 Common uncertainty metric

All portfolios must use the **same physical uncertainty normalization**.

Define a global block structure, e.g.

\[
\Delta_{\rm all}
=
\operatorname{blkdiag}
(\Delta_{30},\ldots,\Delta_{38})
\]

with fixed \(W_i^L,W_i^R\).

Portfolio \(S\) selects the corresponding blocks.

Never fit a different uncertainty scale separately for every portfolio and then compare the resulting radii.

If repair/controller changes \(u\) alter the validity domain, either:

- use a single envelope valid on the entire design domain, or
- use declared parameter-dependent weights \(W(s;\theta,u)\).

## 2.11 Robust transition structure

Define base radius:

\[
r_{\rm base}=r_\Delta^\perp(\emptyset).
\]

Only for:

\[
0\le\varepsilon<r_{\rm base}
\]

define the standard simplicial compatibility complex:

\[
\mathcal K_\varepsilon(\theta)
=
\{
T:
r_\Delta^\perp(R;\theta)>\varepsilon
\ \forall R\subseteq T
\}.
\]

Define minimal robust blockers:

\[
\mathcal H_\varepsilon(\theta)
=
\min_{\subseteq}
\{
S:r_\Delta^\perp(S;\theta)\le\varepsilon
\}.
\]

Define any-order transition margin:

\[
m_\Delta^\perp(T;\theta)
=
\min_{S\subseteq T}
r_\Delta^\perp(S;\theta).
\]

Blocker persistence:

\[
b_H=r_\Delta^\perp(H),\qquad
d_H=\min_{G\subsetneq H}r_\Delta^\perp(G)
\]

including \(G=\emptyset\).

Then:

\[
H\in\mathcal H_\varepsilon
\iff
b_H\le\varepsilon<d_H.
\]

Do not oversell this as deep new combinatorics; it is structural organization.

---

# 3. CANONICAL NUMBERS THAT MUST BE REPRODUCED BEFORE NEW SCIENCE

The following are **Gate 0 references**, not assumptions to force-fit.

## 3.1 V4 / P4

Candidate core:

```text
{30, 33, 35, 37}
```

At frozen P4:

```text
g = 0.03625
k = 1.425
t = 1.5
h = 1
```

Expected full H4 result:

```text
alpha = +0.1270065 s^-1
f = 0.6222797 Hz
```

All 15 proper subsets stable.

Expected H4 controller boundary:

```text
g* ≈ 0.20768
f* ≈ 0.706 Hz
```

Do not confuse P4 interior with the boundary.

## 3.2 Exact closure audit

Expected:

```text
minimum physical local singular value ≈ 0.2973
collective sigma_min(I + Q_H) ≈ 1.26e-8
minimum proper-subset collective factor ≈ 0.2945
min |lambda(Q_H) + 1| ≈ 3.51e-8
max contextual-return distance to +1 ≈ 2.28e-7
Schur residuals < 3.7e-16
```

## 3.3 V9 targeted family

Expected same-policy P4 target-band result:

```text
512 portfolios
511/512 classification
1 false-safe
0 false-unstable
14 minimal blockers
5 blockers order 4
6 blockers order 5
3 blockers order 6
kappa = 4
predicted and full-order minimal blocker antichains identical
```

Expected maximum reduced/full blocker root errors approximately:

```text
real part <= 2.72e-7 s^-1
frequency <= 8.23e-9 Hz
```

## 3.4 V9 full spectrum caveat

Expected:

```text
402/512 correct
110 false-safe
0 false-unstable
all 110 false-safe are aperiodic-real in the archived final audit
```

The targeted EM construction is **not** a full-spectrum safety certificate.

## 3.5 Policy dependence

Expected documented qualitative behaviour includes:

```text
4 -> 3 -> 2 -> 3 -> 4 -> empty
```

on a frozen F7B line.

Expected distinct hypergraph counts reported historically:

```text
30 / 36 / 16
```

across the three policy planes.

## 3.6 Governed model

Expected P4 H4 changes approximately:

```text
+0.1270 -> -0.0745 s^-1
```

with TGOV1N governors.

Do not hide this. It is evidence that the blocker is policy/model conditioned.

## 3.7 TDS

Expected archived phasor-domain evidence includes:

```text
32/32 declared verdicts
80/80 P4 Monte-Carlo verdicts
measured rate vs alpha_perp correlation r ≈ 0.969
```

## 3.8 Condenser / retuning

Expected:

```text
Q/V Newton -> g* ≈ 0.20768 in 5 iterations, ~2e-10 final error
damped condenser threshold ≈ 2.48% = 106.0 MVA
TDS threshold ≈ 2.47%
```

Also preserve the separate finite-disturbance envelope requirement:

```text
~684–792 MVA
```

for the declared 200 MW disturbance envelope. It is **not** the small-signal threshold.

## 3.9 Topology campaign

Known retrospective results:

```text
P4 branch doublings 0, 1, 13, 43 independently remove H4
42 policy/action pairs create an incompatibility of size <= 3
169/810 policy-action pairs change H
```

These are no longer eligible as “blind” results.

## 3.10 Line-sensitivity transfer caveat

Expected historical result:

```text
equation-equivalent model: 12/12 signs, Spearman 1.00
cross-converter transfer: 7/12 signs
```

This caveat is scientifically important and motivates robust/model-holdout work.

---

# 4. REQUIRED DIRECTORY TREE FOR THE NEW CAMPAIGN

Create all new work under one root, e.g.

```text
research/ias2026_bulletproof/
├── README.md
├── docs/
│   ├── PREREGISTRATION.md
│   ├── CLAIM_LEDGER.md
│   ├── NEGATIVE_RESULTS.md
│   ├── DEVIATIONS.md
│   ├── MODEL_SCOPE.md
│   ├── POWERDYNAMICS_RECONCILIATION.md
│   ├── SECOND_MODEL_SPEC.md
│   ├── ROBUST_UNCERTAINTY_SPEC.md
│   ├── REVIEWER_ATTACK_MATRIX.md
│   └── FINAL_EXECUTIVE_REPORT.md
├── prereg/
│   ├── policies.json
│   ├── blind_predictions.json
│   ├── uncertainty_weights_manifest.json
│   ├── topology_holdout.json
│   └── sha256.txt
├── env/
│   ├── python/
│   │   ├── requirements.txt
│   │   └── python_env.txt
│   ├── julia/
│   │   ├── Project.toml
│   │   ├── Manifest.toml
│   │   └── julia_env.txt
│   └── system_info.txt
├── src/
│   ├── python/
│   ├── julia/
│   ├── shared/
│   └── tests/
├── raw/
│   ├── gate0/
│   ├── powerdynamics/
│   ├── ports/
│   ├── robust/
│   ├── tds/
│   ├── topology/
│   ├── holdout/
│   └── scaling/
├── derived/
│   ├── tables/
│   ├── npz/
│   ├── jld2/
│   └── manifests/
├── figures/
│   ├── png/
│   ├── pdf/
│   └── svg/
├── logs/
├── reports/
│   ├── GATE0_REPORT.md
│   ├── PD39_GATE_REPORT.md
│   ├── PORT_RECONCILIATION_REPORT.md
│   ├── ROBUSTNESS_REPORT.md
│   ├── HOLDOUT_REPORT.md
│   ├── SCALING_REPORT.md
│   ├── TDS_REPORT.md
│   ├── ABLATION_REPORT.md
│   └── FINAL_BULLETPROOF_AUDIT.md
└── paper/
    ├── main.tex
    ├── main.pdf
    └── figures/
```

Never write new outputs into historical TX4/CDW result directories.

---

# 5. PHASE 0 — ENVIRONMENT + PROVENANCE

## 5.1 Repository inventory

Produce:

```text
reports/REPOSITORY_INVENTORY.md
```

Include:

- git root;
- HEAD;
- branch;
- relevant tags;
- all known historical TX4/CDW branches;
- status;
- canonical result directories;
- Python version;
- Julia version;
- OS;
- CPU;
- BLAS;
- package versions.

Save:

```bash
git status --short
git log --all --decorate --oneline --max-count=100
git branch -a
git tag --list
```

## 5.2 Python environment

Prefer an isolated environment.

Example:

```bash
python -m venv .venv-ias2026
# activate appropriately for OS
python -m pip install --upgrade pip wheel setuptools
python -m pip install numpy scipy pandas matplotlib sympy networkx \
    pyyaml h5py jupyter pytest statsmodels scikit-learn \
    control cvxpy
python -m pip freeze > research/ias2026_bulletproof/env/python/requirements.txt
```

If the existing repository already has a lockfile/environment, use it first and record deviations.

## 5.3 Julia environment

Do not modify a global Julia environment.

Create:

```text
research/ias2026_bulletproof/env/julia/
```

Activate it.

At minimum investigate/install:

```julia
PowerDynamics
ControlSystemsBase
RobustAndOptimalControl
CairoMakie
DataFrames
CSV
JLD2
LinearAlgebra
SparseArrays
OrdinaryDiffEqRosenbrock
OrdinaryDiffEqNonlinearSolve
Graphs
```

Use:

```julia
using Pkg
Pkg.activate("research/ias2026_bulletproof/env/julia")
Pkg.instantiate()
Pkg.precompile()
Pkg.status()
```

If packages are missing, add them explicitly and commit `Project.toml` + `Manifest.toml`.

Important current capability:

`RobustAndOptimalControl.robstab` currently supports **diagonal complex perturbations**. Do not pretend this is a general arbitrary full-block \(\mu\) implementation. If the intended uncertainty is not representable, use an appropriate conservative alternative and label it correctly.

## 5.4 PowerDynamics sources

Use two sources:

1. current `PowerDynamics.jl` IEEE-39 tutorial (Part I–IV);
2. current `JuliaEnergy/IEEE39.jl` repository if available.

Record exact package/repository commits.

Do not assume an alpha repository is authoritative merely because it is new.

---

# 6. PHASE 1 — GATE 0: REPRODUCE THE FROZEN PYTHON SCIENCE

**Goal:** prove that the current environment can reproduce the canonical nominal results before adding anything.

Create one driver:

```text
src/python/run_gate0.py
```

It must execute or call existing canonical code, never recompute with a new undocumented implementation.

## 6.1 Required outputs

Create:

```text
derived/tables/gate0_v4.csv
derived/tables/gate0_boundary.csv
derived/tables/gate0_v9.csv
derived/tables/gate0_tds.csv
derived/tables/gate0_topology.csv
```

## 6.2 Required pass conditions

Examples:

```text
P4 H4 alpha error <= 1e-6 s^-1
P4 H4 frequency error <= 1e-6 Hz
all 15 proper subset signs identical
g* error <= 1e-5
collective sigma_min order <= 1e-6 at boundary
correct Q sign: -1
correct return sign: +1
V9 target antichain exactly 14 blockers
V9 blocker set equality exact
```

Use reasonable tolerances if historical data precision is lower, but preregister before running.

## 6.3 Failure protocol

If any headline result does not reproduce:

- mark `GATE0 = FAIL`;
- stop any downstream phase that assumes that result;
- compare source code, environment, data files, model parameters;
- do not “update the expected result”.

---

# 7. PHASE 2 — MATHEMATICAL UNIT TESTS / ADVERSARIAL THEORY TESTS

Create:

```text
src/python/tests/test_port_identities.py
src/python/tests/test_contextual_return.py
src/python/tests/test_bridge.py
src/python/tests/test_transverse.py
src/julia/test_theory.jl
```

## 7.1 Random complex-matrix tests

At least 10,000 randomized cases across dimensions 1–8.

Test:

1. determinant lemma;
2. local/collective factorization;
3. Schur contextual return;
4. \(q_i\neq0\) under theorem assumptions;
5. bridge:
   \[
   G=K(I+DK)^{-1}=K(I+Q)^{-1}L^{-1};
   \]
6. robustness bounds;
7. basis changes;
8. near-singular local factors;
9. singular \(D\);
10. nonnormal examples.

Required residual targets:

```text
median relative residual < 1e-13
max residual < 1e-10 unless deliberately ill-conditioned
```

For ill-conditioned tests, report condition number and backward error rather than silently relaxing tolerance.

## 7.2 Counterexample search

Randomly search for counterexamples to every theorem whose assumptions can be violated.

Examples:

- singular \(L_S\);
- improper subset singular;
- repeated critical root;
- mode-band crossing without imaginary-axis crossing;
- ill-conditioned \(A_S,C_S\);
- hidden unobservable unstable mode.

Document the failures as demonstrations of why assumptions are necessary.

---

# 8. PHASE 3 — POWERDYNAMICS IEEE-39 BASELINE GATE

This phase is **independent-code validation**.

Historical warning: a previous PowerDynamics campaign ended in `CASE C — STOP` because apparent blockers were not reproducible under equilibrium-path audits. Therefore this phase begins with the operating point, not blockers.

## 8.1 Build the official/tutorial baseline

Use the current PowerDynamics IEEE-39 tutorials.

Run:

- model creation;
- initialization;
- power flow;
- dynamic initialization;
- small-signal linearization if exposed;
- small-disturbance simulation.

## 8.2 Initialization path audit

Run at least:

```text
tolerances: 1e-8, 1e-10, 1e-12
multiple initial guesses
multiple solver methods if practical
power-flow initialized route
direct dynamic initialization route if available
```

For each solution save:

```text
equilibrium residual
bus V magnitudes
bus angles
generator P/Q
dynamic state values
Jacobian condition estimates
critical eigenvalues
```

## 8.3 Gate A pass condition

The operating point must be path-consistent.

Pre-register a tolerance such as:

```text
max |ΔV| <= 1e-7 pu
max |Δtheta| <= 1e-7 rad
max critical-mode ΔRe <= 1e-4 s^-1
max critical-mode Δf <= 1e-4 Hz
```

If PowerDynamics cannot meet a justified tolerance:

```text
PD39_GATE_A = STOPPED_BY_GATE
```

Do not claim independent blocker validation.

You may still use PowerDynamics for other non-dependent experiments, but explicitly label them.

## 8.4 Baseline modal reconciliation

Compare Python frozen base and PowerDynamics base:

- inter-area mode frequencies;
- damping;
- mode shapes;
- machine participation;
- MAC-like mode similarity after mapping states.

Do not require identical spectra if the device models differ.

Produce:

```text
figures/.../pd39_base_spectrum.*
figures/.../pd39_mode_mac.*
derived/tables/pd39_base_modes.csv
```

---

# 9. PHASE 4 — SAME-MODEL CROSS-CODE VALIDATION

**Purpose:** prove the TX4 phenomenon is not a Python implementation artifact.

Implement the frozen GFL equations in Julia as literally as possible.

## 9.1 Component-by-component unit tests

Before placing it in IEEE-39, build a one-bus/infinite-bus harness.

Test independently:

- PLL;
- P/Q filter;
- outer P/Q PI;
- inner current PI;
- Q/V loop;
- terminal current injection;
- per-unit scaling.

At an identical equilibrium, compare Python vs Julia:

```text
state derivative
finite-difference Jacobian
analytic/AD Jacobian
terminal admittance over frequency
```

Frequency grid:

```text
0.01 Hz – 100 Hz
log spaced >= 400 points
dense around 0.2–2 Hz
```

Target tolerances for exact same-model port:

```text
relative transfer-matrix error median < 1e-6
max < 1e-4 excluding numerically singular points
```

If conventions differ (dq orientation, current sign, Park transform), fix the mapping explicitly and document it.

## 9.2 Same-model PD39 V4

Install the same GFL at candidate buses 30,33,35,37 using matched dispatch and ratings.

Run all 16 portfolios.

Compare:

- equilibrium;
- \(\alpha_\perp\);
- target mode;
- Hasse lattice;
- blockers.

There are two possible scientific outcomes:

### Outcome A
H4 or an equivalent blocker reproduces closely.

Label:

```text
POWERDYNAMICS_VALIDATED
```

### Outcome B
The exact blocker differs but a collective blocker exists.

Still valuable. Continue mechanism validation.

### Outcome C
No reproducible blocker exists.

Record negative result. Do not force it.

## 9.3 Same-model mechanism

For any reproducible PowerDynamics blocker:

- build exact terminal port;
- compute local factors;
- compute \(Q_H\);
- continue to boundary;
- verify \(-1\);
- verify contextual return \(+1\);
- compare port root to full DAE root.

This is more important than reproducing the exact set `{30,33,35,37}`.

---

# 10. PHASE 5 — ALTERNATIVE CONVERTER MODEL HOLDOUT

This phase tests whether the mechanism survives a materially different converter model.

## 10.1 Do not invent a fake “second model”

Priority:

1. documented model already available in PowerDynamics or another public Julia implementation;
2. published phasor GFL model with equations and parameters;
3. generic model based on a documented standard.

A GFM droop model is **not** a substitute for a GFL holdout. It can be a separate technology-extension experiment.

## 10.2 Match equilibrium, not dynamics

For each SG replacement:

- match scheduled \(P,Q\);
- match MVA rating;
- solve the same network equilibrium;
- use the alternative controller dynamics unchanged except required operating setpoints.

## 10.3 Discovery then blind holdout

Split policies deterministically before evaluating full-order labels.

Example:

```text
DISCOVERY:
5 policy points

HOLDOUT:
8 policy points generated and saved before any full eigensolve
```

The reduced/port construction may inspect terminal models and predict blockers on holdout.

Freeze predictions to:

```text
prereg/blind_predictions.json
```

with SHA-256 and git commit.

Only after freeze run full DAE.

Metrics:

```text
exact blocker-antichain match
blocker precision
blocker recall
minimum-order kappa match
root Re error
root frequency error
false-safe count
false-unstable count
```

This is the genuine blind experiment.

---

# 11. PHASE 6 — NONLINEAR TDS VALIDATION IN JULIA

Use common small disturbances designed for small-signal consistency.

Do not use a severe fault as the sole validation of a small-signal theorem.

## 11.1 Required cases

For at least:

```text
base
one stable proper subset
one blocker
same blocker after repair/retune
one governed case
one alternative-model case
```

## 11.2 Disturbance

Use one common small disturbance across compared cases, e.g.:

- +2% bus-20 active-load pulse;
- or another frozen disturbance from the original campaign.

Keep equal horizons.

## 11.3 Signals

Collect:

- relative machine speeds;
- GFL estimated frequencies;
- selected angle differences;
- bus voltage magnitudes;
- currents.

## 11.4 Estimators

Implement at least two:

1. matrix pencil / Prony;
2. bandpass + Hilbert/log-envelope fit.

Require agreement within preregistered tolerance.

Compare:

\[
\alpha_{\rm eig}
\]

vs measured envelope rate.

Figures:

```text
eig_vs_tds_scatter
common_disturbance_traces
mode_frequency_comparison
```

---

# 12. PHASE 7 — TRANSVERSE ROBUST LFT AND PHYSICAL UNCERTAINTY ENVELOPE

This phase is required before claiming robust transition compatibility on IEEE-39.

## 12.1 Construct the transverse nominal system

Use either:

```text
A_perp route
```

or:

```text
symmetry-deflated exact port route
```

Verify no structural zero remains.

## 12.2 Obtain model discrepancy

Preferred discrepancy source:

### Route A — same operating point, alternative GFL
For every candidate bus:

\[
E_i(j\omega)
=
Y_i^{alt}(j\omega)
-
Y_i^{nom}(j\omega).
\]

### Route B — EMT vs phasor
If a trustworthy EMT terminal scan exists, use it.

### Route C — parameter ensemble
Only if A/B unavailable. Label as parameter uncertainty, not EMT uncertainty.

## 12.3 Build a GLOBAL envelope

Fit \(W_i(s)\) once over the declared domain.

Required coverage plots:

```text
singular values of E_i(jw)
weight magnitude
normalized residual
coverage over all calibration models/parameters
```

Define the physical meaning of:

```text
epsilon = 1
```

e.g. “one full calibrated terminal-model envelope”.

Do not use arbitrary percent language unless physically justified.

## 12.4 RobustAndOptimalControl usage

If using `robstab` / `structured_singular_value`, document that the current package implementation supports diagonal complex perturbations.

Only claim exact \(\mu\) for a block structure actually supported.

For unsupported full blocks:

- use fixed-frequency unstructured singular-value radius;
- use conservative diagonal embedding;
- or implement/validate a separate upper-bound calculation.

Label every number.

---

# 13. PHASE 8 — REAL V4 ROBUST CENSUS

This is the decisive new experiment.

For all 16 real V4 portfolios compute:

```text
nominal alpha_perp
nominal target mode
r_delta_perp(S)
worst-case frequency
active uncertainty block
sigma_min(I+Q_S) at relevant frequency
local factor conditioning
bridge lower/upper bounds
```

Output:

```text
derived/tables/v4_robust_radii.csv
```

## 13.1 Core robust objects

Compute:

\[
\mathcal H_\varepsilon^\perp,
\quad
I_H,
\quad
\Pi_H,
\quad
m_\Delta^\perp(T).
\]

Use epsilon grid selected before inspecting the final atlas.

## 13.2 Key falsification question

Does the robust radius collapse when the collective closure approaches singularity?

Do not merely plot correlation.

For each stable portfolio / policy point calculate:

```text
r_delta_perp
sigma_min(I+Q)
bridge lower bound
bridge upper bound
worst-case frequency
closure-minimum frequency
```

Test:

- whether the worst-case robust frequency matches or tracks the collective-closure frequency;
- when the bounds become loose;
- whether looseness is explained by \(A_S,C_S\) conditioning.

This is a much stronger validation of the Collective–Robustness Bridge than random matrices.

---

# 14. PHASE 9 — REAL CONTROLLER × UNCERTAINTY ATLAS

Use a real frozen controller path, not a synthetic policy coordinate.

Recommended first path:

```text
P4 Q/V gain g
```

including the known boundary near:

```text
g* ≈ 0.20768
```

Choose a preregistered g-grid.

For every \(g\):

1. solve equilibrium;
2. evaluate all 16 V4 portfolios;
3. compute \(r_\Delta^\perp\);
4. construct \(\mathcal H_\varepsilon^\perp(g)\);
5. compute \(m_\Delta^\perp(T;g)\).

Required flagship plot:

```text
(g, epsilon) -> minimum robust blocker order
```

Second plot:

```text
surface/heatmap of m_delta_perp(T;g)
with color or annotations for active bottleneck coalition
```

Third:

```text
blocker persistence barcode at selected g values
```

This is the figure that determines whether “Robust Transition Compatibility” becomes a headline or remains future work.

---

# 15. PHASE 10 — TOPOLOGY + CONTROL REPAIR

Separate retrospective validation from blind prediction.

## 15.1 Retrospective topology calibration

Known repairs:

```text
branch doublings 0,1,13,43
```

Use them only for:

- mechanism analysis;
- calibration;
- comparing metrics.

Never label them blind.

For all 46 branches compare:

```text
Fiedler
Kirchhoff
gSCR
minSCR
direct eigenvalue sensitivity
total re-equilibrated sensitivity
closure sensitivity
robust-transition margin change
finite full-DAE outcome
```

## 15.2 Total derivative

For physical topology change \(\rho\):

\[
F(\xi,\rho)=0
\]

\[
\xi_\rho=-F_\xi^{-1}F_\rho
\]

and:

\[
\widehat{\mathcal T}_\rho
=
\mathcal T_\rho
-
D_\xi\mathcal T[
F_\xi^{-1}F_\rho].
\]

Define reproducibly:

```text
direct structural
power-flow/network/load-mediated
device/controller-mediated
```

via an explicit partition of \(\mathcal T\).

## 15.3 New blind topology holdout

Create a genuinely unseen task.

Good choices:

- governed holdout policy;
- alternative converter holdout;
- unseen reinforcement magnitude;
- unseen topology action set.

Before reveal save predictions:

```text
top-k stabilizing actions
sign of effect
predicted blocker removed/created
predicted change in m_delta_perp
```

Then reveal via full DAE + TDS where appropriate.

## 15.4 Any-order repair

Optimization objective should be physical or explicitly a surrogate.

If using a finite action library, say:

> lowest-cost feasible action **within the fixed candidate library**

not global minimum.

Compare:

### Final-only design
satisfy final portfolio margin.

### Any-order design
satisfy:

\[
m_\Delta^\perp(T;u)\ge\varepsilon_{\rm req}.
\]

Use the same physical cost definition.

Report:

```text
final margin
worst intermediate margin
active bottleneck coalition
cost
new blockers created
```

---

# 16. PHASE 11 — GOVERNOR HOLDOUT

Our P4 blocker disappears with documented governors. This is not an embarrassment; it is a test of mechanism generality.

Search the governed policy plane for a **reproducible minimal blocker**.

Do not choose it after looking for a visually pleasing result.

Predefine selection rule, e.g.:

```text
lowest-cardinality blocker with minimality separation >= 0.02 s^-1
then largest separation
then lexicographic tie-break
```

For the chosen governed blocker:

- verify full DAE;
- TDS;
- local regularity;
- \(Q_H\to-1\);
- return \(+1\);
- one repair.

If no governed blocker passes:

```text
GOVERNED_HOLDOUT = NEGATIVE
```

and report it.

---

# 17. PHASE 12 — ABLATION CAMPAIGN

Run ablations to answer “what part is actually necessary?”

Each ablation must have a question and outcome metric.

Required ablations:

1. **No local dressing**
   \[
   Q \to DK_o
   \]
   Does prediction degrade?

2. **Static network kernel**
   freeze \(K(j\omega)\) or use low-frequency/static approximation.

3. **No controller self-energy**
   use simplified L/M/D-like model.

4. **No Q/V loop**

5. **No PSS**

6. **With vs without governor**

7. **Matched MW/MVA controls**

8. **Frozen equilibrium vs re-equilibrated topology**

9. **Single-device sensitivity vs pairwise vs exact portfolio**

10. **Target band vs complete spectrum**

11. **Same-model vs alternative-model GFL**

12. **Nominal repair vs robust any-order repair**

For each produce:

```text
exact hypergraph match
precision/recall
root error
ranking Spearman
sign accuracy
false-safe rate
```

---

# 18. PHASE 13 — SCALABILITY BENCHMARK

Do not use a larger bus system as proof of physical generality unless dynamic data are documented.

## 18.1 IEEE-118 network

Julia can load MATPOWER `case118`.

Use it for scalability if no documented dynamic case is available.

Explicit label:

```text
SCALABILITY BENCHMARK WITH PARAMETERIZED DYNAMICS
```

not independent physical validation.

## 18.2 Candidate-set sizes

Run:

```text
m = 6, 8, 10, 12
```

and if feasible:

```text
m = 15
```

Compare:

- full enumeration;
- port reduced evaluation;
- active-set / blocker-cut search;
- memory;
- wall time.

Report:

```text
number of portfolios evaluated
number of blockers discovered
runtime
speedup
peak memory
```

Do not claim polynomial scaling if enumeration remains exponential.

## 18.3 Larger documented dynamic case

If a documented 145/162/etc. dynamic benchmark can be imported cleanly, use it as an additional holdout.

Provenance matters more than bus count.

---

# 19. PHASE 14 — FIGURES REQUIRED FOR THE FINAL PAPER / POSTER

Every figure must be generated programmatically.

Export:

```text
PNG >= 300 dpi
vector PDF
SVG
```

Use consistent fonts and labels.

## Core figures

### F1 — V4 Hasse lattice
Show all 16 portfolios, stable proper subsets, H4 unstable.

### F2 — P4 vs boundary
Small table/diagram distinguishing:

```text
P4:
g=0.03625, alpha=+0.1270, f=0.6223 Hz

boundary:
g*=0.20768, alpha=0, f≈0.706 Hz
```

### F3 — local vs collective singular values
Show:

```text
local min ≈ 0.2973
collective -> 1.26e-8
proper collective >= 0.2945
```

### F4 — complex-plane closure
Plot critical eigenvalue of \(Q_H\) approaching \(-1\) and contextual returns approaching \(+1\).

### F5 — policy geometry
Show blocker/hypergraph transitions, not just kappa.

### F6 — V9 blocker antichain
Display all 14 target-family blockers, by order and membership.

### F7 — TDS
Proper subset decays, blocker grows, repaired blocker decays.

### F8 — PowerDynamics cross-code
Mode/eigenvalue correspondence and MAC/participation.

### F9 — uncertainty-envelope calibration
Nominal vs alternative/EMT terminal maps + weights.

### F10 — real V4 robust radii
All 16 real portfolios.

### F11 — real \((g,\varepsilon)\) robust compatibility atlas
Potential award figure.

### F12 — Collective–Robustness Bridge
Scatter/trajectory:

```text
sigma_min(I+Q)
vs
r_delta_perp
```

plus theoretical bounds.

### F13 — topology repair map on IEEE-39
Color branches by physically meaningful action effect.

### F14 — final-only vs any-order repair
Show margin of every intermediate portfolio and physical cost.

### F15 — blind holdout
Predicted vs revealed blocker/action results.

### F16 — scalability
Runtime vs candidate count / portfolio count.

---

# 20. PHASE 15 — TABLES REQUIRED

Generate at minimum:

```text
TABLE_MODEL_PARAMETERS.csv
TABLE_V4_PORTFOLIOS.csv
TABLE_V9_BLOCKERS.csv
TABLE_POLICY_HYPERGRAPHS.csv
TABLE_CLOSURE_AUDIT.csv
TABLE_TDS.csv
TABLE_POWERDYNAMICS_RECONCILIATION.csv
TABLE_SECOND_MODEL_HOLDOUT.csv
TABLE_UNCERTAINTY_WEIGHTS.csv
TABLE_V4_ROBUST_RADII.csv
TABLE_REPAIR_ACTIONS.csv
TABLE_ABLATIONS.csv
TABLE_SCALING.csv
TABLE_CLAIM_LEDGER.csv
```

Each table must contain source-file provenance.

---

# 21. PHASE 16 — STATISTICAL / NUMERICAL HYGIENE

## 21.1 Never report only a correlation

For ranking experiments report:

- Spearman;
- sign accuracy;
- top-k recall;
- precision-recall if positives are rare;
- rank of known effective actions;
- uncertainty / bootstrap CI where meaningful.

## 21.2 Conditioning

For every critical root report:

\[
|\ell^H\mathcal T_s r|
\]

and a clearly defined normalization-dependent condition metric.

Do not call it coordinate invariant.

## 21.3 Finite differences

For every derivative validation run a step-size sweep.

Example:

```text
1e-2, 1e-3, 1e-4, 1e-5, 1e-6
```

Show a convergence plateau.

## 21.4 Root tracking

Use continuity / modal assignment, not nearest-frequency only.

Track with a combination of:

- eigenvalue distance;
- mode-shape MAC;
- participation.

Log ambiguous crossings.

---

# 22. PHASE 17 — REVIEWER ATTACK MATRIX

Create:

```text
docs/REVIEWER_ATTACK_MATRIX.md
```

For every central claim include:

```text
claim
strongest possible objection
evidence
remaining limitation
status
what would falsify it
```

Mandatory attacks:

1. artifact of omitted governor;
2. artifact of custom GFL;
3. artifact of Python implementation;
4. artifact of operating-point path;
5. artifact of selected mode band;
6. failure outside 0.3–1.5 Hz;
7. EMT mismatch;
8. structural zero breaks robust LFT;
9. uncertainty metric not comparable across portfolios;
10. repair exits calibration domain;
11. sensitivity not new;
12. minimal failure sets not new;
13. small benchmark;
14. exponential scaling;
15. current limits / DC link absent;
16. topology rankings fail model transfer;
17. no blind evidence;
18. no physical cost model.

Do not write marketing responses. Write technically precise answers.

---

# 23. PHASE 18 — CLAIM LEDGER

Create CSV and Markdown versions.

Columns:

```text
claim_id
claim_text
mathematical_status
evidence_status
benchmark
scope
source_files
novelty_status
reviewer_risk
allowed_poster_wording
forbidden_wording
```

Example:

```text
C01
At the frozen P4 benchmark, all 15 proper V4 portfolios are stable while H4 is unstable.
IEEE39_VALIDATED
...
Allowed:
"At the frozen IEEE-39 P4 benchmark..."
Forbidden:
"Four GFLs are generally incompatible."
```

No claim enters the final abstract unless:

```text
reviewer_risk <= medium
```

or the limitation is stated in the same sentence.

---

# 24. PHASE 19 — PAPER REWRITE RULES

Only after all experiments are frozen.

## 24.1 Abstract

Lead with real IEEE-39 evidence.

Do **not** lead with toy robust numbers.

Required order:

1. problem;
2. MRB;
3. exact collective closure;
4. IEEE-39 flagship;
5. V9 antichain / policy dependence;
6. independent validation / robust result if it actually passes;
7. one sentence of scope.

## 24.2 Results hierarchy

Use:

```text
Tier 1: exact/proved theory
Tier 2: frozen IEEE-39 results
Tier 3: independent PowerDynamics/model holdout
Tier 4: robust physical uncertainty results
Tier 5: synthetic explanatory pilots
```

Never reverse this hierarchy.

## 24.3 Limitations

Explicitly state:

- phasor-domain;
- no current limits;
- no DC-link dynamics;
- no modulation delay;
- no ride-through logic;
- exact robust structure supported/not supported by Julia package;
- second-model status;
- EMT status.

---

# 25. PHASE 20 — POSTER DESIGN SPEC

The poster must be understandable at three distances.

## 25.1 Three meters

One headline:

> EVERY SMALLER TRANSITION IS STABLE. THE COMPLETE TRANSITION IS NOT.

One V4 lattice.

## 25.2 One meter

Four blocks:

```text
DISCOVER -> EXPLAIN -> GENERALIZE -> ACT
```

### DISCOVER
MRB / Hasse.

### EXPLAIN
local regular vs collective closure.

### GENERALIZE
policy geometry + V9.

### ACT
retuning + condenser + topology + robust transition design if validated.

## 25.3 Thirty centimeters

Only here show:

- exact equations;
- contextual return;
- robustness bridge;
- limitations / QR reproducibility.

---

# 26. PHASE 21 — FINAL ZIP FORMAT

At the end create exactly one archive:

```text
IAS2026_BULLETPROOF_CLOSURE_<YYYYMMDD_HHMM>.zip
```

Root structure:

```text
IAS2026_BULLETPROOF_CLOSURE/
├── README_FIRST.md
├── FINAL_EXECUTIVE_REPORT.pdf
├── FINAL_EXECUTIVE_REPORT.md
├── FINAL_PAPER.pdf
├── FINAL_POSTER_DRAFT.pdf
├── CLAIM_LEDGER.csv
├── CLAIM_LEDGER.md
├── NEGATIVE_RESULTS.md
├── OPEN_LIMITATIONS.md
├── REPRODUCIBILITY_MANIFEST.json
├── SHA256SUMS.txt
├── git/
│   ├── git_status.txt
│   ├── git_log.txt
│   ├── branch_and_head.txt
│   └── commits_created.txt
├── prereg/
├── env/
├── src/
│   ├── python/
│   └── julia/
├── data/
│   ├── raw/
│   └── derived/
├── figures/
│   ├── png/
│   ├── pdf/
│   └── svg/
├── reports/
├── paper_source/
└── logs/
```

## 26.1 ZIP size

If raw traces are huge:

- use `.csv.gz`, `.npz`, `.jld2`;
- include all summary-level data;
- include a manifest for any omitted huge trace;
- prefer a ZIP under a practical upload size.

## 26.2 README_FIRST.md

Must contain:

1. scientific verdict;
2. which gates passed/failed;
3. headline claims safe to use;
4. claims that failed;
5. exact commands to reproduce;
6. key file paths;
7. environment versions;
8. final branch/HEAD;
9. no-push statement if applicable.

---

# 27. REQUIRED EXECUTION COMMANDS / ORCHESTRATORS

Create top-level runners.

## Python

```text
src/python/run_all_python.py
```

CLI examples:

```bash
python src/python/run_all_python.py --phase gate0
python src/python/run_all_python.py --phase theory-tests
python src/python/run_all_python.py --phase robust-v4
python src/python/run_all_python.py --phase topology
python src/python/run_all_python.py --phase figures
python src/python/run_all_python.py --phase all
```

## Julia

```text
src/julia/run_all_julia.jl
```

CLI examples:

```bash
julia --project=research/ias2026_bulletproof/env/julia \
  src/julia/run_all_julia.jl gate-pd39

julia --project=research/ias2026_bulletproof/env/julia \
  src/julia/run_all_julia.jl same-model

julia --project=research/ias2026_bulletproof/env/julia \
  src/julia/run_all_julia.jl robust-v4

julia --project=research/ias2026_bulletproof/env/julia \
  src/julia/run_all_julia.jl tds
```

Every phase must be restartable and checkpointed.

No phase should require rerunning hours of unrelated computation after interruption.

---

# 28. AUTOMATIC PASS/FAIL DASHBOARD

Generate:

```text
reports/GATE_DASHBOARD.csv
reports/GATE_DASHBOARD.md
```

Rows:

```text
G0 Frozen Python reproduction
G1 Theory unit tests
G2 PowerDynamics equilibrium
G3 Same-model Julia/Python ports
G4 Same-model blocker mechanism
G5 Alternative-model holdout
G6 Julia TDS
G7 Transverse robust LFT
G8 Common uncertainty envelope
G9 Real V4 robust census
G10 Real policy×uncertainty atlas
G11 Retrospective topology calibration
G12 Genuine blind topology/model holdout
G13 Governed blocker
G14 Ablations
G15 Scalability
G16 Final claim audit
```

Columns:

```text
status
metric
threshold
observed
evidence_file
notes
```

The final report must show this dashboard on page 1.

---

# 29. EXPECTED DECISION LOGIC

## If PowerDynamics reproduces same-model blocker

Excellent. Promote cross-code evidence.

## If PowerDynamics has a different blocker but same collective mechanism

Also excellent. Promote mechanism generality, not set equality.

## If PowerDynamics has no blocker under physically stable initialization

Do not hide it. State that blocker existence is model/policy conditioned; use PowerDynamics as a negative holdout.

## If alternative GFL changes blocker

Expected and valuable.

## If alternative GFL removes all blockers

Still valuable; strengthens model-dependence caveat.

## If robust radius does not track collective closure

The Collective–Robustness Bridge theorem may still be algebraically correct but not practically predictive because the prefactors/uncertainty directions dominate. Report this.

## If robust atlas is trivial

Do not promote “Robust Transition Compatibility” as headline.

## If robust atlas is rich and physically calibrated

Promote it.

---

# 30. FINAL SCIENTIFIC SUCCESS CRITERIA

The work is “bulletproof enough for the poster” only if all of the following are true:

### A. Nominal phenomenon
At least one frozen, reproducible minimal blocker exists in the canonical model.

### B. Exact mechanism
For at least one blocker:

```text
local factors regular
Q_H -> -1
contextual return -> +1
```

with numerical residuals commensurate with conditioning.

### C. Nonlinear consistency
Small-disturbance TDS agrees in sign/frequency/rate within declared tolerance.

### D. Independent evidence
At least one of:

```text
PowerDynamics same-model
alternative converter
documented dynamic benchmark
```

supports either the blocker or the mechanism.

### E. Blind evidence
At least one genuinely unseen prediction is frozen and revealed.

### F. Robustness
If robustness is a headline:
a physically interpretable common uncertainty metric must be used on real V4 portfolios.

### G. Actionability
At least one control/support/topology intervention has a finite full-order validation.

### H. Scope honesty
Full-spectrum, EMT, current-limit, DC-link, and model-transfer limitations are explicit.

If D/E/F fail, the nominal blocker poster may still be excellent, but robust-generalization claims must be reduced.

---

# 31. FINAL INSTRUCTIONS TO SONNET

Do not stop after producing plots.

Your final responsibility is to answer, with evidence:

1. Which original claims survived?
2. Which claims were weakened?
3. Which were refuted?
4. Does PowerDynamics reproduce the phenomenon or only the mechanism?
5. Does a second converter model reproduce the mechanism?
6. Is the robust radius physically meaningful and comparable across portfolios?
7. Does collective closure predict robustness collapse in the real IEEE-39 case?
8. Can a repair improve the entire transition without creating a smaller blocker?
9. Does the method predict a genuinely unseen case?
10. What is safe to put on the Vancouver poster?
11. What must remain supplementary/future work?
12. What is the strongest, most defensible novelty statement after all negative results?

Then generate the final ZIP.

Do **not** claim completion until:

```text
SHA256SUMS.txt exists
all figure source data exist
all plots are regenerated from code
all main numbers appear in a source CSV/JSON
the claim ledger is complete
the negative-results file is complete
the final ZIP opens and lists correctly
```

At the very end print exactly:

```text
FINAL BUNDLE:
<absolute path to zip>

FINAL BRANCH:
<branch>

FINAL HEAD:
<commit>

GATES PASSED:
<list>

GATES FAILED/UNRESOLVED:
<list>

TOP 5 SAFE POSTER CLAIMS:
1.
2.
3.
4.
5.

TOP 5 FORBIDDEN/UNSUPPORTED CLAIMS:
1.
2.
3.
4.
5.
```

---

# 32. DO NOT DO THESE THINGS

- Do not use known branches 0,1,13,43 as a “blind” reveal.
- Do not call random synthetic matrices IEEE-39 validation.
- Do not call closure-space \(r_Q\) a physical terminal uncertainty radius.
- Do not use a \(\mu\) routine outside its supported uncertainty structure and then call it exact.
- Do not compare robustness radii computed under different uncertainty normalizations.
- Do not say “PowerDynamics validates” if its equilibrium gate fails.
- Do not call a GFM holdout a second GFL model.
- Do not claim EMT validation unless an EMT gate actually passes.
- Do not hide the 110 full-spectrum V9 aperiodic false-safe cases.
- Do not hide that governors stabilize the P4 flagship.
- Do not call synthetic L1 action magnitude “engineering cost”.
- Do not call the 106 MVA condenser threshold a finite-disturbance support requirement.
- Do not treat the physical switching transient between configurations as analyzed if only stationary re-equilibrated configurations were studied.
- Do not optimize the story after seeing holdout labels.
- Do not delete failed results.

---

# 33. PRIORITY ORDER IF COMPUTE/TIME BECOMES LIMITED

If resources are constrained, execute in this order:

```text
P0  Gate 0 frozen reproduction
P1  Theory/transverse corrections
P2  PowerDynamics equilibrium gate
P3  Same-model Julia/Python port reconciliation
P4  Same-model PowerDynamics V4 + mechanism
P5  Julia TDS
P6  Alternative-model blind holdout
P7  Common uncertainty envelope
P8  Real V4 robust census
P9  Real g×epsilon atlas
P10 Physical any-order repair
P11 Ablations
P12 IEEE-118/large-case scalability
```

Never sacrifice P0–P6 merely to obtain a large-bus figure.

---

# 34. ONE-SENTENCE TARGET

The final campaign should allow the paper/poster to say something close to the following **only if the data support it**:

> A generation transition is a family of dynamically distinct configurations: inclusion-minimal replacement portfolios can fail through an exact network-mediated collective return while every smaller portfolio remains regular; this mechanism can be identified from terminal dynamics, tested across independent implementations and converter models, and used to design controller, synchronous-support, or topology actions that protect the entire implementation path rather than one final eigenvalue.

If the independent/robust gates do not support the full sentence, shorten it rather than stretching the evidence.
