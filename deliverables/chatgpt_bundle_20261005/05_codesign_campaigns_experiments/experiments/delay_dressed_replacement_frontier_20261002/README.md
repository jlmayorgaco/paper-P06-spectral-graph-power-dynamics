# Delay-dressed graph damping and SG→GFL replacement frontier

This directory contains a new, isolated experiment against the repository's IEEE-39 model. It preserves the frozen historical experiment outputs and records the source state before this directory was created. No commit or push was made.

**Decision: the requested delayed maximum-replacement result is not established.** The no-delay parity gate passed. The exact retarded characteristic matrix, local root continuation, local sensitivities, frozen delay patterns, and lossless Taylor validation are present. Two independent gates block the proposed result: complete DDE rightmost-root coverage was not established, and the gauge-deflated hidden block for the proposed exact Schur damping operator is singular at s=0. Therefore there is no validated D_G decomposition, no delay co-design, and no delay-dependent maximum or global/near-global optimality claim.

Start with [STATUS.md](STATUS.md), [MODEL_PROVENANCE.md](MODEL_PROVENANCE.md), and [DELAY_MODEL_CONTRACT.md](DELAY_MODEL_CONTRACT.md). The numerical derivation is in [THEORY_DELAY_REPLACEMENT.md](THEORY_DELAY_REPLACEMENT.md); actual supported statements only are in [POSTER_CLAIMS.md](POSTER_CLAIMS.md). Reproduction steps and known block are in [REPRODUCE.md](REPRODUCE.md).

Key outputs are under `baseline_reproduction/`; pattern vectors and graph operators are at the experiment root. Files in `inputs/` are immutable snapshots of prior parameter inputs, not newly reproduced data. DDE results contain no Padé approximation.

The complete folder is bundled as `delay_dressed_replacement_frontier_20261002.zip`. Its archive includes the source, inputs, frozen vectors, diagnostics, reports, tables, figures, and SHA-256 inventory.
