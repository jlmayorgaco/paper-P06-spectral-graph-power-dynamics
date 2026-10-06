# TX4 Gate-0 Reference Manifest

**Canonical uploaded artifact:** `TX4_FINAL_MODAL_SCOPE_CHATGPT_UPLOAD.zip`  
**SHA-256:** `548b259c51264917c8bb13deddc633dde74e884ef0f4f73667d83e8b8dc48982`  
**Archive entries:** 33  
**Role:** authoritative *reference-output bundle* for the final TX4 modal-scope audit.  
**Not sufficient by itself for reproduction:** the archive contains docs, result tables and figures, but not the complete executable TX4 source tree / environment.

## Provenance declared by the bundle

- Branch: `research/tx4-final-modal-scope-audit`
- Starting commit: `2cbb860eb4ba053edee99c1a040c8024f3f6828e`
- Parent freeze: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Content-artifact commit: `1380413d`
- Push: none
- Final case: `A`

## Gate-0 authoritative numerical targets

### V4 / P4

- Portfolio: `30+33+35+37`
- Full-order transverse alpha: `0.127006467144 s^-1`
- Frequency: `0.622279669578 Hz`
- All 15 proper V4 subsets: stable in `TX4_V4_BLIND_VS_FULL.csv`
- Only V4 nonnegative-alpha row: `30+33+35+37`

### Correct collective-return sign audit

- nearest `Q_H` distance to `-1`: `3.509875453894e-08`
- maximum contextual-return distance to `+1`: `2.280537928238e-07`
- maximum Schur residual: `3.625891133257e-16`

Correct convention:

```text
Q_H -> -1
R_i|R -> +1
```

### Physical local vs collective audit

From `TX4_TRUE_LOCAL_VS_COLLECTIVE.csv`:

- minimum physical `sigma_min(I+M_ii)` over the stored sweep:
  `0.297326880959346`
- minimum stored collective `sigma_min(I+Q_H)`:
  `1.261484530209848e-08`

Headline reference additionally reports:

- minimum proper-subset collective sigma at the boundary frequency:
  `0.294546050542390`

### V9 same-policy P4 scope

From `TX4_V9_GLOBAL_VS_EM_TRUTH.csv`:

- global false-safe: `110`
- EM-band false-safe: `1`
- EM-band false-unstable: `0`
- aperiodic global modes: `271`
- global modes in target band: `233`

Authoritative headline counts:

```text
global: 402/512 correct, 110 false-safe, 0 false-unstable
target 0.3–1.5 Hz: 511/512 correct, 1 false-safe, 0 false-unstable
```

The 110 global false-safe cases are aperiodic-real misses.

### TDS traces in this upload

Three equal-horizon / archived comparison cases are present:

| case | eig alpha | eig f | fitted/primary alpha | primary f | verdict |
|---|---:|---:|---:|---:|---|
| proper 30+33+35 | -0.204688 | 0.916019 | -0.204998 | 0.916076 | STABLE |
| H4 original 30+33+35+37 | 0.127006 | 0.622280 | 0.127358 | 0.622232 | UNSTABLE |
| H4 retuned g=0.25 | -0.017385 | 0.710541 | -0.017414 | 0.710546 | STABLE |

## Important traps / non-authoritative fields

### 1. Do not reuse normalized local columns in `TX4_CONTEXTUAL_RETURN_BOUNDARY.csv`

That file contains:

```text
local_sigma_min_at_device = 1
local_factor_det_abs = 1
```

for all four device rows. Those are not the corrected physical local factors `I+M_ii`.

The corrected physical-local authority is:

```text
results/TX4_TRUE_LOCAL_VS_COLLECTIVE.csv
```

and the final report explicitly states that normalized `I+Q` diagonal blocks are identity blocks and cannot establish physical local regularity.

### 2. Two slightly different collective-sigma values appear in the bundle

`TX4_CONTEXTUAL_RETURN_BOUNDARY.csv` stores approximately:

```text
2.241536e-08
```

while the final physical sweep / headline stores:

```text
1.261484530209848e-08
```

Treat the latter as the authoritative *minimum over the corrected frozen sweep*, and do not silently equate the two.  
If a new campaign needs the exact boundary-point value, regenerate it from the executable canonical source and document why the two archived evaluations differ (e.g. root/evaluation interpolation or different fresh audit path).

### 3. `TX4_FINAL_MODAL_SCOPE_DEVIATIONS.md` contains duplicated `Deviation D1`

This is a packaging/editorial duplication, not two independent deviations. Preserve it as historical evidence; do not count it twice.

### 4. Legacy global V9 numbers are explicitly non-authoritative for same-policy P4

Do not use:

```text
395/512
116 false-safe
```

as same-policy P4 evidence. The bundle states that those came from an FC10/default path.

### 5. This archive is evidence, not a runnable repository

The archive contains only:

```text
docs/
figures/
results/
```

It does not contain the full model source, environment lockfiles, experiment scripts, or raw state-space construction needed to reproduce Gate 0 from scratch.

## Files to treat as canonical authorities

### Sign / return
- `results/TX4_Q_VS_RETURN_SIGN_AUDIT.csv`
- `docs/TX4_CONTEXTUAL_RETURN_THEOREM.md`

### Physical local-vs-collective
- `results/TX4_TRUE_LOCAL_VS_COLLECTIVE.csv`
- `results/TX4_TRUE_COLLECTIVE_SUBSET_FACTORS.csv`

### V4 truth
- `results/TX4_V4_BLIND_VS_FULL.csv`

### V9 same-policy modal/global truth
- `results/TX4_V9_GLOBAL_VS_EM_TRUTH.csv`
- `results/TX4_V9_EM_BAND_CONFUSION.csv`
- `results/TX4_V9_FALSE_SAFE_TAXONOMY.csv`
- `results/TX4_V9_SLOW_FALSE_SAFE_AUDIT.csv`

### TDS included in this bundle
- `results/TX4_TDS_FINAL.csv`
- `results/TX4_TDS_SUMMARY.json`

### Claim wording
- `results/TX4_FINAL_MODAL_SCOPE_CLAIM_MATRIX.csv`
- `docs/TX4_FINAL_MODAL_SCOPE_REPORT.md`
- `docs/TX4_FINAL_MODAL_SCOPE_CHATGPT_HANDOFF.md`

## What Sonnet still needs for the full bulletproof campaign

One of the following is required:

1. the actual repository at the canonical freeze / relevant later branch; or
2. a reproducibility ZIP that includes the model source, experiment scripts, dependencies and raw inputs.

Minimum missing source categories:

```text
TX4 DAE/model implementation
SG/GFL model equations in executable form
V4/V9 portfolio runners
port-construction code
boundary-continuation code
TDS runner
policy-map scripts
topology-action scripts
environment/requirements
```

Without those, Sonnet can audit the evidence bundle but cannot honestly claim to have *reproduced* it.

## Gate-0 instruction

When the executable repository is available:

1. compute new outputs in a new directory;
2. never read these CSVs as computational inputs;
3. compare the new outputs against this manifest only after execution;
4. preserve all mismatches;
5. fail Gate 0 if a headline result cannot be reconciled within preregistered tolerance.

