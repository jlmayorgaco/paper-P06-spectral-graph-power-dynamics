# P6 genuine holdout v2

The input file `p6_genuine_holdout_inputs.csv` freezes four discovery gains and four unseen holdout gains for the official PowerDynamics GFL11 alternative-model harness. Discovery labels are generated first. Predictions for the four holdout gains are then written and hashed without reading holdout labels. Only after that hash is recorded is the holdout phase executed.

This protocol is an alternative-model blind holdout, not the true same-model custom SynchronousMachine/AVR/PSS gate. The holdout is valid only for the declared PowerDynamics model and common H4 portfolio. If the revealed holdout is single-class, balanced accuracy, unstable precision/recall, MCC, and Cohen kappa are reported as not meaningful and a second preregistered holdout is required before any discrimination claim.
