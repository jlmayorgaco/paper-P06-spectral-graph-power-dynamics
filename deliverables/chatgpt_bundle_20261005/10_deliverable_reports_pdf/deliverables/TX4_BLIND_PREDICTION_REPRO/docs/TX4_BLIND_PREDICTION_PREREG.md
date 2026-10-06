# TX4 blind portfolio prediction and contextual-return preregistration

Version: `TX4-BP-v1`  
Frozen before prediction/evaluation: yes  
Parent: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`  
Branch: `research/tx4-blind-portfolio-prediction-final`

## Scope and stop rules

The scientific question is whether one reusable baseline network kernel and
local SG-to-GFL replacement models can predict the minimal incompatible TX4
portfolio before consulting the full-order portfolio answer table, and whether
the exact contextual return explains its clean g-only boundary. Existing
full-order outputs remain read-only and are loaded only after immutable blind
prediction artifacts are committed.

The campaign has three statuses: PASS, FAIL, and BLOCKED. A blocked descriptor,
ANDES boundary, or EMT test is never turned into a positive claim. No push is
authorized.

Forbidden claims and campaigns: Schur complement, determinant, Nyquist,
return-ratio, eigenvalue sensitivity, or Newton as new mathematics; Shapley,
cumulants, cycle holonomy, Fiedler/Laplacian mechanism, generic weakness,
hidden radius, generic co-design, EMT portfolio validation, universal
generalization, or “first” without a focused literature audit.

## Frozen model and answer-key firewall

V4 candidates are `{30,33,35,37}`; V9 candidates are the existing frozen
`{30,31,32,33,34,35,36,37,38}` set. The policy is P4
`(g,k,t,h)=(0.03625,1.425,1.5,1)`. The clean path varies only `g` with
`(k,t,h)=(1.425,1.5,1)`. The search band is `0.3--1.5 Hz` positive
frequency.

Before the reveal, prediction code may read only: the frozen model equations;
the all-SG baseline; the network Ybus and operating-point data; local SG/GFL
device definitions; controller parameters; candidate bus/rating definitions;
and the exact port-closure implementation. It may not query any file or
function containing full-order portfolio alpha/verdict labels, H, kappa,
historical boundary values, F4/V4 tables, or V9 tables. Prediction outputs are
written to new files, hashed, and committed. The answer-key comparison occurs
in a separate post-reveal command.

## Blind V4 protocol

Enumerate all 16 V4 subsets. For each subset, build `D_S`, `M_S`, and `Q_S`
from the one baseline kernel and local Delta-Y library. Search the fixed
positive-frequency band for closure roots using the frozen frequency grid and
root refinement: 121 uniform frequencies in `0.3--1.5 Hz`, then a fixed
bounded local refinement around the best candidates. Use no full-order
portfolio eigenanalysis. Label a portfolio `UNSTABLE_PREDICTED` only when a
right-half-plane closure root is found with regular local factors; otherwise
`STABLE_PREDICTED` or `UNRESOLVED_PREDICTED`. Save root, frequency, closure
residual, local-factor minimum, and runtime.

The primary V4 gate is false-safe count zero; the preferred gate is 16/16
verdict agreement after reveal. Prediction artifacts must be committed before
the reveal.

## Blind V9 protocol

Using the same baseline K and local Delta-Y library, enumerate all 512 V9
subsets and write immutable predictions before opening the V9 full-order
answer table. Use the same root protocol and no answer-key screening. After the
prediction commit, compare against the frozen full-order V9 table. Report
false-safe/false-unstable counts, exact H and kappa recovery only where the
answer table defines them, and unresolved rows.

## Full-order reveal and runtime

After each blind prediction commit, load the historical full-order answer key
and report verdict accuracy, alpha/root error, frequency error, H, and kappa.
Full-order eigenanalysis is the ground truth, not a competitor. Measure full
and reduced runtimes with `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`
and the same hardware/runtime. The reduced measurement includes one-time
baseline/local setup and per-portfolio closure evaluation; the full
measurement includes direct full model/equilibrium/Jacobian/transverse
eigenanalysis. If reduced is not faster, report that.

## Contextual return and minimality

With `C_H=I+Q_H`, `R=H\\{i}`, and

```
G_R=(I+Q_RR)^(-1)
R_i|R=Q_iR G_R Q_Ri,
det(I+Q_H)=det(I+Q_RR)det(I-R_i|R).
```

The theorem is the classical block determinant identity under regularity,
nonsingular local factors, minimality, nonsingular proper-subset closure, and
a simple collective root. At the full/reduced root report every return
eigenvalue nearest 1, `|mu-1|`, `sigma_min(I-R)`, `cond(I+Q_RR)`,
`sigma_min(I+Q_RR)`, local factors, `eta_H`, `sigma_H`, and `tau_H`.

Sweep the fixed g grid and locate the return unity minimum without using the
full-order g boundary. Agreement target is `|g_return-g_full|<=1e-4` and
`max|mu-1|<=1e-6` or a declared solver floor. Local factors must remain
regular while the collective factor closes.

## Aggregate, screens, and remediation

Before viewing alpha differences, select same-cardinality V9 pairs by exact
MW/MVA match if available, otherwise nearest normalized MW/MVA distance with
the fixed lexicographic tie-break `(distance, subset labels)`. A stable/unstable
pair with aggregate mismatch <=5% is a supporting result; otherwise report
NONE.

Reuse frozen singleton/pair/triple, modal, additive, pairwise, third-order
alpha, determinant-truncation, aggregate, and static results. Keep third-order
alpha distinct from third-order determinant truncation. Compare them with full
DAE and exact reduced closure. Do not say full eigenanalysis fails.

For remediation use only `g`; report the reduced direction, first-order
prediction, root refinement, and fixed stable point `g_after=0.25`. Do not
claim economic optimality.

## TDS, ANDES, and EMT

Use the frozen G2 IEEE-39 D2 active-load pulse at bus 20, +2%, duration 0.2 s,
with the archived BDF/guards/observables. Run H4 at P4, H4 at `g_after`, and
the four P4 triples if runtime permits. This is nonlinear phasor-domain TDS,
not EMT. Reuse the existing ANDES validation without retuning; an unsupported
g-boundary comparison is `NOT TESTED`. EMT portfolio validation is
`UNRESOLVED`, with no optional stretch run in the core campaign.

## Literature, claims, and case rule

Audit 2020--2026 primary literature on grid strength/gSCR, impedance and
gain-phase stability, stability manifolds/operating-point sets, root-cause,
placement, coordination, robust stability, minimal cut sets, and discrete
replacement portfolios. Avoid “first/unique/universal” unless directly
supported.

Case A requires H4 minimality, 15 proper subsets, local regularity, blind V4
prediction, material reduced frequency/root agreement, contextual return,
g-boundary agreement, and TDS consistency; V9/runtime/aggregate matches are
bonuses. Case B means the narrow flagship story survives but prediction,
runtime, aggregate, or independence is limited. Case C means the flagship or
collective mechanism fails or is blocked.
