# Model contract — frozen IEEE-39 TX4/H4 mechanism campaign

Status: `PHASE0_INVENTORY_COMPLETE`

This document freezes the model and numerical gates before any new campaign
code is written. The campaign must import the existing implementation and may
add only routing, measurement, and reporting adapters under this campaign
directory.

## 1. Scientific target

- Network: frozen IEEE-39 case, 100 MVA, 60 Hz, constant-power loads.
- Replaced converter buses: `H4 = (30, 33, 35, 37)`.
- Frozen parameter point: `P4 = (g, k, t, h) = (0.03625, 1.425, 1.5, 1.0)`.
- GFL policy: the frozen leaky voltage-control law used by `_f7_common.py`,
  with `LEAK = 0.05`.
- Requested mechanism: exact DAE/index-1 reduction, transverse modal
  reduction, retained port operator, normalized closure operator, exact Schur
  complement/contextual-return factors, local-versus-collective separation,
  and root-motion sensitivity.

The target gate values are those already registered by the TX4 modal-scope
campaign:

| gate | expected value |
|---|---:|
| H4 state dimension | 86 dynamic states |
| H4 transverse rightmost target alpha | `+0.127006468 s^-1` |
| H4 target frequency | `0.6222796695 Hz` |
| boundary gain | `g* ~= 0.20768` |
| proper-subset status at P4 | all stable |
| H4 status at P4 | unstable |
| physical local-factor minimum at P4 | `0.2973268809593455` or larger |
| collective boundary singular value | `~1.26e-8` |

The model contract treats the sign convention as fixed: `det(I+Q_H)=0`
corresponds to the `-1` eigenvalue of `Q_H`, while a contextual return reaches
`+1`. The physical local factor is the pre-normalized diagonal factor
`I+M_ii`, not the diagonal of a normalized `I+Q` matrix.

The local-factor gate is evaluated at the frozen P4 point and at the
eigenvalue boundary. A wider exploratory gain grid is also recorded, but a
lower local factor away from the registered P4/boundary mechanism points is
reported as an observation rather than silently folded into the canonical
minimality claim.

## 2. Reusable frozen implementation

The active Python implementation is imported from the committed package below;
it is not copied into the run directory:

```text
reports/poster/ias2026/research/src/ibr_cycles/
reports/poster/ias2026/research/experiments/_f7_common.py
reports/poster/ias2026/research/experiments/tx4_contextual_return.py
```

Relevant entry points:

- `ibr_cycles.models.ieee39_case.ReplacementPlan`, `solve_case`:
  equilibrium solve, central-difference Jacobians, and index-1 reduction.
- `ibr_cycles.models.ieee39_devices.ConverterParameters`:
  frozen converter parameterization.
- `ibr_cycles.models.port_admittance.build_action_space`:
  retained port/action operator.
- `ibr_cycles.certification.symmetry.rotation_generator` and
  `frequency_partner`: rotational/transverse mode construction.
- `ibr_cycles.certification.transverse.transverse_operator`:
  quotient transverse operator.
- `_f7_common.solve_subset`, `Theta`, `LEAK`: registered TX4 parameter path.

The exact reusable model files currently have these SHA-256 hashes:

```text
models/ieee39_case.py       1f38a8e541780ab7205809dc1ba58bdafb6c923933f5772d46407064edc8c722
models/ieee39_devices.py    bcf0c76be65b171b580684065cde4eb88d0433c6b4d7beefec966f8711036d90
models/ieee39_network.py    f42144cd02a6dc779ed690fd866ad026f33eb7844584c7f10281fc7d98e3fe13
models/port_admittance.py   b0fa694dd0111be158133a9cf8c4e36be670e50ef562dcd63f89d7ff87b8d4ef
models/port_core.py         38c54b1a70f87fbbe2c301a717e6dd6b4112b3a9cf8deb663426ce08b729bb9b
experiments/tx4_contextual_return.py
                              77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74
```

The frozen network source config is
`research/ias2026_last_validation/raw/true_same_model/canonical_source/configs/ieee39_network.json`
with SHA-256
`705f53299b57cbd536ababfe5236f3b2ec621e8b305a5737135cfba8363d995e`.
The referenced `data/raw/ieee39_full.xlsx` has SHA-256
`9c2048dc94201ee48ffe816de65d53831fa1f4b0e7fdf32f24c50f1017edcfe5`.

## 3. Independent Julia path

The independent same-math Julia reference environment is retained at:

```text
research/ias2026_last_validation/code/julia/run_true_same_model_ieee39.jl
research/ias2026_last_validation/code/julia/frozen_gfl11.jl
research/ias2026_last_validation/env/julia/Project.toml
research/ias2026_last_validation/env/julia/Manifest.toml
```

The available environment was verified as Julia `1.11.9`, PowerDynamics
`5.0.0`, NetworkDynamics `1.3.0`, ModelingToolkitBase `1.74.0`, SciCompDSL
`1.0.4`, OrdinaryDiffEqNonlinearSolve `2.9.8`, SciMLBase `3.55.0`,
ControlSystemsBase `1.22.0`, DataFrames `1.8.2`, Graphs `1.15.0`, and
CairoMakie `0.15.14`. Julia is an independent parity/reference path; it must
not silently replace the frozen Python TX4 model.

The historical `P2_POWERDYNAMICS_V4_STATUS.md` and
`P4_JULIA_TDS_STATUS.md` are alternative-model/TDS evidence, not canonical
same-model proof.

## 4. Historical evidence and model-selection warning

Registered evidence is kept under:

```text
research/ias2026_last_validation/raw/modal_scope/
research/ias2026_last_validation/reports/
research/ias2026_last_validation/raw/true_same_model/
```

The modal-scope evidence reports the required H4 value
`alpha = +0.1270064671440848 s^-1` and frequency
`0.6222796695779256 Hz`. A separate older `raw/true_same_model` census reports
an H4 transverse value near `0.14467 s^-1`; that is retained as historical
provenance but is explicitly a distinct snapshot/convention for this campaign.
It must not be mixed into the target gate. If the active source cannot
reproduce the registered `0.127006468` gate with the frozen P4 definition,
the campaign is `BLOCKED` and no tuning or model substitution is allowed.

## 5. Repository and run hygiene

- Git branch policy: only local `main`; no push is authorized.
- Inventory HEAD before campaign work: `d931adb4`.
- New artifacts belong only under
  `reports/poster/ias2026/research/bnd_h4_mechanism/results/<RUN_ID>/`.
- Each run uses an immutable run id
  `YYYYMMDDTHHMMSS_<git-shortsha>_<tag>` and contains its own `raw/`,
  `derived/`, `tables/`, `figures/`, `claims/`, `report/`, `environment/`,
  and `logs/` directories.
- No new `.py`, `.jl`, `.csv`, `.png`, `.md`, `.zip`, or logs may be written
  to the repository root or loose campaign siblings.
- Historical scripts that hard-code `ROOT / results` must not be executed
  directly. The campaign will use a routing adapter or an isolated subprocess
  whose output root is the immutable run directory.
- Existing historical files are preserved; this campaign does not delete or
  overwrite them. No archive/zip is produced for this run.

## 6. Reproducibility gate

Before reporting science, the run manifest must record the exact commit,
working-tree state, Python/Julia versions, source hashes, command lines, and
all output paths. The campaign must report gates G0--G11 explicitly. Any
failed canonical gate is reported as `BLOCKED` with the observed values and
the run is not promoted to a scientific claim.
