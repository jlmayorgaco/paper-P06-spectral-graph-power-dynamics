# BND / IEEE-39 research bundle — read this first (2026-10-05)

Repository: paper-P06-spectral-graph-power-dynamics (branch `codex/collective-interaction-bounds-20261003`, HEAD 8f21ec0b, plus uncommitted work).
Author: Jorge Luis Mayorga Taborda, Universidad de los Andes. Target: IEEE IAS 2026 student poster competition (Vancouver).
Topic: replacing synchronous generators (SG) with inverters (GFL, PLL-based) in IEEE-39 and IEEE-68 — "Beyond Nodal Damping" (BND).

**Start with `01_STATE_OF_RESEARCH.md`, then `02_CLAIM_STATUS_AND_NOVELTY.md`, then `03_PROMPT_FOR_CHATGPT.md`.**
`MANIFEST.csv` lists every file with SHA-256; `SKIPPED_FILES.csv` lists what was left out (too large / binary).

## Folder map
| Folder | What is inside |
|---|---|
| `01_papers_series` | LaTeX/PDF of the paper series: TX1 controller-aware bridge, TX2 "when is weak a graph object", TX3 connected intervention calculus, TX4 policy-dependent incompatibility, IAS minimal dynamic incompatibility, CDW contextual dynamic weakness. |
| `02_portfolios_incompatibility_TX3_TX4` | Portfolio / incompatibility line: `spectral_portfolio_theory_v1`, PD39 portfolio campaign results (255+1 closure, blocker atlas, confirmatory, robust transition), TX3 artifacts, TX4 final reports and claim matrices. |
| `03_weak_nodes_CDW_algebraic` | Weak nodes / weak links with algebraic characterisation: Contextual Dynamic Weakness (CDW) campaign, hardening, IEEE-68 replication (Case B), TX2 weak-as-graph-object paper, weak-elements conference draft. |
| `04_BND_H4_mechanism_and_research_line` | H4 mechanism study (portfolio H4={30,33,35,37} instability), research docs/theory/results of the IAS-2026 research tree. |
| `05_codesign_campaigns_experiments` | All `experiments/*` of Oct 1-5 2026: analytic iteration, nonlinear co-design, delay-dressed frontier, exact action space (rank<=30), mega analytical delay co-design, latency-robust PLL, collective damping theory, interaction decision (certified counterexample + repair), network-moment laws, graph gain geometry, frequency limits / physical bridge, material synthesis, theory synthesis, predictive frequency bounds (protocol only). |
| `06_codesign_campaigns_reports_AtoQ` | Earlier experiment series A-Q (`reports/experiment_*`): spectral co-design, Jordan/mixed, secure GFL KKT optimum (Q2B, 94.2% candidate later retracted). |
| `07_poster_sources_and_outputs` | Poster LaTeX sources (versions 20261003 and 20261004), generated figures, compiled PDFs/PNGs (`output/pdf`). |
| `08_handoffs_and_state` | README, architecture, reproducibility, runbooks, previous ChatGPT handoff documents, AGENTS.md. |
| `09_model_source_code` | Julia model source (`src/pd39`), Project.toml, scripts. |
| `10_deliverable_reports_pdf` | Compiled reports in `deliverables/`. |

## Deliberately NOT included
Large traces (>120 KB CSV), `.npz/.zip/.jld2` binaries, build directories, Julia depots. Hashes in `MANIFEST.csv` identify what is here.
The directory `series_replacement_portfolios/` (Parts I-VI of a six-part series referenced in the author's notes) was NOT found in this repository; the series papers present are those under `01_papers_series`.
The term "LD/LDM" in the request was not identifiable; everything under the repo's result directories is included by category instead.
