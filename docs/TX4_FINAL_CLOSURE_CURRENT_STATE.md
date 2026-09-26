# TX4 final-closure current state

This file records the archaeological state before the contextual-return
campaign. It is not a new numerical result.

## Frozen parent and branch

- Parent commit: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Parent tag: `TX4_FINAL_MANUSCRIPT_FREEZE`
- Campaign branch: `research/tx4-contextual-return-final`
- Clean worktree: `C:\\w\\tx4closure`
- No push is authorized.

The parent is the final TX4 manuscript freeze. The dirty working tree in the
user's main checkout is not used by this campaign.

## Frozen TX4 witness

The authoritative implementation is the IEEE-39 full-order phasor-domain DAE
under `reports/poster/ias2026/research/src/ibr_cycles/`. The canonical point is

```
P4 = (g=0.03625, k=1.425, t=1.5, h=1)
H4 = {30, 33, 35, 37}
```

The frozen FC01 structure file reports `alpha_perp(H4)=+0.1270064666836382`
s^-1, all 15 proper subsets with negative transverse spectral abscissa, a
two-dimensional transverse quotient, and two right-half-plane transverse
modes. The frozen FC01 spectral identity residual is approximately `1.55e-8`.
H4 is not described as an irreducible four-device interaction.

The clean controller-only path fixes `(k,t,h)=(1.425,1.5,1)` and varies only
`g`. The archived PCV04/FC13 result places the H4 crossing at

```
g*=0.20768140450381395
alpha approximately 0
frequency=0.7064247879345361 Hz
d alpha/d g=-0.46100603930620704 s^-1
```

These are provenance anchors. The new campaign must independently re-evaluate
them before using them as validation outcomes.

## Exact code and data paths

- Model/equilibrium: `reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py`
- Device equations: `reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py`
- Network: `reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_network.py`
- Exact port action space: `reports/poster/ias2026/research/src/ibr_cycles/models/port_admittance.py`
- Direct transverse path: `reports/poster/ias2026/research/experiments/_f7_common.py`
- Frozen nonlinear TDS: `reports/poster/ias2026/research/src/ibr_cycles/nonlinear/tds.py`
- Frozen TDS configuration: `reports/poster/ias2026/research/configs/ias2026/final_nonlinear_composability_v1.yaml`
- Frozen structure audit: `reports/poster/ias2026/research/results/FINAL_CLOSURE/FC01_structure.csv`
- Frozen controller boundary: `reports/poster/ias2026/research/results/PCV/PCV04/PCV04_summary.json`
- Frozen ANDES validation: `reports/poster/ias2026/research/outputs/ias2026/final_validation_overnight_20260910T003225/E31_ANDES_validation/`

## Scope boundary

The campaign narrows TX4 to exact collective closure, contextual return,
g-only boundary continuation, local-versus-collective separation, proper-subset
minimality, g-only remediation, return derivative verification, one common
nonlinear phasor-TDS disturbance before/after, reuse of the existing ANDES
result, a compact conventional-screen comparison, and a focused literature
gap. Weak-node/link ranking, structured radius, repair optimization, planners,
PowerDynamics/PD39, IEEE-68, EMT, and new controller searches are out of scope.

The existing ANDES result is an independent phasor computation with different
SG, exciter, stabilizer, and network assembly. It does not validate custom-GFL
EMT behavior. ParaEMT portfolio validation remains unresolved.
