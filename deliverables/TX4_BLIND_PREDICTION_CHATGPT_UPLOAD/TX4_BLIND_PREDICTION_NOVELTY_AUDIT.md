# TX4 focused literature novelty audit

Date: 2026-09-14. Scope: 2020--2026 primary research and classical methods
relevant to grid strength, impedance/return-ratio analysis, converter
stability regions, placement, and minimal-cut reasoning. This is a focused
positioning audit, not a systematic review. It does not justify the words
“first”, “unique”, “unprecedented”, or “universal”.

| reference | object studied | continuous/discrete decision | one operating configuration / portfolio family | minimal incompatible device support? | exact return characterization? | local controller remediation? | overlap | distinction |
|---|---|---|---|---|---|---|---|---|
| Cao et al., IEEE TSG 2020, [DOI](https://doi.org/10.1109/TSG.2020.2978250) | Low-frequency multi-inverter islanded microgrids | Continuous terminal characteristics and frequency-domain analysis | Multi-inverter systems | No finite SG-to-GFL replacement minimality test identified | Yes, return-ratio / return-difference GNC | Not the TX4 finite-support test | Multiport return analysis | TX4 asks whether a finite set of local replacements is the minimal unstable support in a full IEEE-39 DAE and conditions one device on the remaining set |
| Zhang et al., IEEE TIE 2022, [DOI](https://doi.org/10.1109/TIE.2021.3095791) | Grid-connected inverter PLL/grid-impedance interaction | Controller and impedance shaping | Single/multi-inverter weak-grid studies | No | Impedance ratio / Nyquist | Yes | Grid-following synchronization interaction | TX4’s H4 is a discrete portfolio witness with a contextual return extracted from the same network kernel |
| Wang et al., IEEE TPWRS 2024, [DOI](https://doi.org/10.1109/TPWRS.2023.3319708) | Slow converter-driven stability with GFL/GFM devices | Continuous parameter bounds and system studies | Two- and five-inverter examples | No finite minimal replacement support reported | Return-ratio/modal analysis | Yes | Multi-device interactions and system-level planning | TX4’s contribution is the finite subset question plus blind V4 prediction and minimality separation |
| Coordinating IBRs, IEEE TSTE 2025, [DOI](https://doi.org/10.1109/TSTE.2024.3406758) | Small-signal stability-constrained IBR operation/coordination | Continuous active-power coordination | Test-system operating scenarios | No | Grid-strength-oriented small-signal constraints | Yes | Planning under stability constraints | TX4 does not claim a superior optimizer; it isolates a minimal discrete incompatibility and tests reduced closure prediction |
| Conte et al., “Small-Signal Stability Manifolds”, arXiv 2026, [record](https://arxiv.org/abs/2605.26254) | Stability regions/manifolds in converter-dominated systems | Continuous controller/operating parameters | Many operating scenarios | Not the finite replacement-support question | Full-network eigenanalysis and learned boundary approximation | Yes | Parameter-space stability boundaries | TX4 uses an exact port factorization for a fixed model and a discrete portfolio family; the g-boundary is a model-specific validation, not a new stability-manifold theorem |
| Caro-Ruiz et al., SEGAN 2020, [DOI](https://doi.org/10.1016/j.segan.2020.100302) | Minimum cut-set vulnerability in power networks | Discrete topological attacks | Cascading / DC-flow network studies | Topological cut sets, not dynamic SG-to-GFL support | No contextual dynamic return | No | Minimal-set language | TX4 uses “minimal” only for a dynamic spectral portfolio property, not reliability cut-set equivalence |

## Answers to the focused questions

The audit found substantial prior work on multi-inverter impedance and
return-ratio analysis, controller-parameter stability regions, and stability-
constrained placement/operation. It did not identify, in the sources checked,
the exact combination of:

1. exhaustive finite SG-to-GFL replacement portfolios;
2. a full-order dynamic failure with every proper replacement subset stable;
3. an immutable pre-reveal reduced-closure prediction of that support; and
4. an exact contextual return conditioned on the remaining devices, validated
   at the same controller boundary.

That is a scoped novelty distinction, not a proof that no prior paper exists.
The mathematical ingredients remain classical: block determinant/Schur
factorization, multiport return ratios, eigenvalue computation, numerical root
finding, and nonlinear phasor-domain TDS. The defensible TX4 novelty is the
experimentally demonstrated packaging and validation of those ingredients for
this finite SG-to-GFL portfolio question, with the V9 negative result kept in
the record.

## Claim discipline

Allowed: “TX4 demonstrates a model-specific minimal dynamic incompatibility
and a reusable reduced closure that predicted the four-bus flagship before
reveal.”

Not allowed: “the first/unique/universal minimal incompatibility,” “a new
Schur theorem,” “EMT validation,” or “generalized stability for all grids.”
