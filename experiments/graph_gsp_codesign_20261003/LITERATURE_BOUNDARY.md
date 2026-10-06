# Literature and novelty boundary

Primary source checked on 2026-10-03:
Linbin Huang, Huanhai Xin, Wei Dong and Florian Doerfler,
*Impacts of Grid Structure on PLL-Synchronization Stability of
Converter-Integrated Power Systems*, arXiv:1903.05489v2 (2019).
https://arxiv.org/abs/1903.05489

That work already uses network Kron reduction, a modal representation, grid
eigenvalue sensitivities, PLL retuning and a 39-bus validation. Graph modes plus
PLL tuning on IEEE-39 are therefore not, by themselves, our novelty claim.
This is a verified nearby comparator, not an exhaustive literature search.

Here the executed question is narrower: with the actual heterogeneous SG/GFL
device dynamics, a frozen 40-ms pure measurement delay and the same local joint
replacement/gain trust box, which gain subspaces preserve all-event feasibility?
The truth model retains the lossy network and complete device dynamics; graph
coordinates restrict the design variables, not the validation model.

Schur complements, Neumann series, commutator identities, implicit root
derivatives and convex LP solutions are established mathematics. We distinguish
their exact application to this model from an original general theorem. A
comparative computational result needs its own evidence. No claim of first-ever
GSP control, superiority to the literature, or global replacement optimality is
made by this pilot.
