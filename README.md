# Spectral Graph Power Dynamics

Research code and editorial artifacts for spectral graph methods in power-system
dynamics.

This repository is organized as a clean architecture research project:

- `src/spectral_ibr/domain`: pure theory objects and services.
- `src/spectral_ibr/application`: use cases and ports.
- `src/spectral_ibr/infrastructure`: adapters for ANDES, solvers, persistence,
  and plotting.
- `src/spectral_ibr/interface`: CLI, configuration, and dependency composition.
- `src/spectral_ibr/experiments`: reproducible gates and manifests.
- `configs`: declarative case, contour, and gate configuration.
- `data`: raw, interim, and processed datasets.
- `outputs`: generated experiment artifacts.
- `reports`: posters, conference papers, Transactions papers, thesis material,
  shared bibliography, and shared notation.
- `temp`: retained scratch area and legacy working material.

`temp/` is intentionally preserved. New reusable code should move into `src/`,
and regenerated scientific artifacts should be written to `outputs/` with a
manifest.

Poster figures should come from gates, not notebooks. The Stage 1 gate builds
processed IBR cases from a citable conversion protocol; poster gates then
consume those processed cases and write versioned artifacts.

## Quick Start

```powershell
python -m pip install -e ".[dev]"
spectral-ibr gate list
spectral-ibr gate run stage1-build-ibr-case --profile mix60
```

Stage 1 reads `configs/cases/conversion_protocol.yaml`. The `base_case` path is
strict: if `data/raw/ieee39_full.xlsx` is missing, the gate fails with an error
instead of searching for a fallback. Change `active_profile` in the config, or
pass `--profile mix20|mix40|mix60|mix80`, to generate:

```text
data/processed/<profile>/<profile>.xlsx
outputs/stage1_build_ibr_case/<profile>/manifest.json
outputs/stage1_build_ibr_case/<profile>/conversion_plan.json
```

After generation, the processed workbook can be loaded by ANDES from Python, for
example:

```python
import andes

system = andes.load("data/processed/mix60/mix60.xlsx", setup=True)
system.PFlow.run()
system.EIG.run()
```

## Research handoff (2026-10-02)

For a continuity package on the IEEE-39 SG-to-GFL PLL co-design and Beyond Nodal Damping work, start with [the next-AI restart guide](reports/handoff/ai_migration_20261002/START_HERE_FOR_NEXT_AI.md), then use the [master research dossier](reports/handoff/ai_migration_20261002/HANDOFF_MASTER.md) and its [artifact map](reports/handoff/ai_migration_20261002/ARTIFACT_MAP.json). The rendered internal report is in output/pdf/IEEE39_BND_PROJECT_HANDOFF_20261002.pdf. This package preserves failed experiments and claim limits; it is not a camera-ready IEEE submission.
