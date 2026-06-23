# arXiv corpus screening memo

Date: 2026-06-23

## Purpose

Build a reproducible arXiv corpus for the thesis literature review on spectral
graph theory, dynamic hosting capacity, rational self-energy, and
converter-dominated power-system dynamics.

This is not yet a final bibliography. It is a large screening corpus that must
be filtered, verified, and connected to peer-reviewed sources before citations
are added to the dissertation or a Transactions-style paper.

## Harvest summary

- API window: 2010-2026
- Topics: 9
- Topic-year query files: 153
- Topic-year records before deduplication: 37,909
- Deduplicated records: 31,477
- Initial ranked PDF sample: 45 PDFs, five per topic
- PDF download errors: 0

Primary data files:

- `data/arxiv_literature/manifest.json`
- `data/arxiv_literature/index/all_metadata.jsonl`
- `data/arxiv_literature/index/all_metadata.csv`
- `data/arxiv_literature/index/summary_by_topic_year.csv`
- `data/arxiv_literature/analysis/screening_report.md`
- `data/arxiv_literature/analysis/top_candidates_by_topic.csv`
- `data/arxiv_literature/pdfs/top_ranked/`

## Topic coverage

| Topic | Records |
|---|---:|
| `krein_sign_dissipation_matrix_polynomials` | 10,069 |
| `spectral_graph_power_networks` | 8,055 |
| `blackbox_converter_modal_identification` | 6,886 |
| `schur_nep_self_energy` | 3,954 |
| `networked_control_resolvent_passivity` | 3,856 |
| `braess_reinforcement_power_networks` | 2,307 |
| `dynamic_hosting_capacity` | 826 |
| `power_system_ibr_stability` | 770 |
| `virtual_inertia_damping_placement` | 186 |

## Immediate reading bins

1. Converter-dominated stability and grid-forming/grid-following dynamics:
   small-signal stability, PLL interaction, weak-grid conditions, GFM/GFL
   placement, and low-inertia behavior.
2. Schur, rational eigenvalue, and self-energy machinery:
   NEPs, transfer-function eigenvalue problems, Schur complements,
   matrix-polynomial sensitivity, and resolvent conditioning.
3. Spectral graph and power-network robustness:
   Laplacian spectra, effective resistance, Kirchhoff index, algebraic
   connectivity, and synchronization robustness.
4. Dynamic hosting and planning:
   hosting capacity, DER hosting, curtailment, protection exposure, and dynamic
   constraints beyond voltage/thermal limits.
5. Reinforcement reversals:
   Braess-like effects, topology control, transmission expansion, line
   reinforcement, and dynamic performance degradation.
6. Control-system certificates:
   passivity, small-gain/small-phase criteria, nonnormality, pseudospectra, and
   decentralized stability.
7. Krein/signature foundations:
   sign characteristics, Krein signatures, dissipation-induced instability,
   indefinite damping, and matrix-polynomial structure.
8. Black-box converter identification:
   impedance models, Loewner/data-driven methods, modal identification, and
   black-box small-signal models.

## Screening rules for the next pass

Each candidate paper should receive:

- `INCLUDE_CORE`: directly supports the thesis argument or paper contribution.
- `INCLUDE_BACKGROUND`: useful context but not central.
- `METHOD_ONLY`: useful method but not a power-system contribution.
- `RELATED_NO_CITE`: adjacent but too far from the thesis claim.
- `EXCLUDE_NOISE`: false positive from broad query terms.

Each included paper should be tagged with:

- model layer: `fixed-D`, `QEP`, `Schur-NEP`, `full-DAE`, `input-output`,
  `graph-only`
- contribution type: `theorem`, `algorithm`, `case-study`, `review`,
  `benchmark`, `negative-result`
- relevance to the main claim: `static-vs-controller sign reversal`,
  `self-energy attribution`, `hosting limit`, `virtual inertia placement`,
  `line reinforcement`, `passivity certificate`, `Krein/signature bridge`

## Next literature-review deliverable

The next writing pass should not simply add more citations. It should rebuild
Chapter 2 around a stronger map:

1. Classical power-system modal stability gives the baseline but not the
   controller-condensed rational operator.
2. IBR/GFM/GFL literature shows the mechanism space where static graph metrics
   can fail.
3. Spectral graph theory provides useful ranking tools but cannot certify
   damping without damping and controller terms.
4. NEP, Schur, and rational eigenvalue literature supplies the mathematical
   language for `Pi(s)`.
5. Passivity and small-gain/small-phase literature supplies complementary
   controller-aware certificates.
6. Krein/sign-characteristic literature is a foundation to adapt, not a
   ready-made theorem for dissipative rational power-system NEPs.
7. The novelty target is the bridge: line/action planning plus Schur
   condensation plus pole-residue attribution plus memory in `T_s`.
