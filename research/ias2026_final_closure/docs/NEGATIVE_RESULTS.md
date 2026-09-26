# IAS2026 negative-results register

This file prevents negative or incomplete evidence from being silently turned
into a positive Vancouver claim.

1. The frozen custom no-governor IEEE-39 model was not ported exactly to Julia.
   The same-model cross-code gate is `STOPPED_BY_GATE`.
2. The PowerDynamics IEEE-39 portfolio census used official machines,
   governors, and AVR devices. Its 16/16 stable result is a valid
   `NEGATIVE_HOLDOUT` for that alternative model, not a refutation of the
   frozen custom model.
3. The SimpleGFLDC IEEE-39 census also produced 16/16 stable portfolios. It is
   another model-conditioned `NEGATIVE_HOLDOUT`, not same-model parity.
4. The official-model TDS cases solved, but the selected trace did not provide
   a reliable oscillation crossing. No modal-frequency match is claimed.
5. The genuine alternative-model holdout was single-class: 4/4 stable.
   Accuracy is reported; balanced accuracy, MCC, kappa, and discrimination are
   undefined. A second preregistered holdout is required before a predictive
   claim.
6. The collective mechanism audit was stopped after the same-model gate.
7. No common physical uncertainty envelope, robust radius, or new scaling
   sweep was executed.
8. No paper or poster was rewritten from these incomplete gates.

These are results, not missing documentation. The raw files and reports keep
the evidence auditable.
