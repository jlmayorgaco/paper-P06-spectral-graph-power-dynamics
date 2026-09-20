# TX4 Exact P4 / GFL11 Julia Reproduction - Current State

Date: 2026-09-19 (America/Bogota)

## Scope

This worktree executes only the frozen TX4 P4/GFL11 cross-language reproduction requested in the attached specification. It does not start a new research campaign, tune parameters, search for new blockers, run SimpleGFLDC, run EMT, or open a theory branch.

## Archaeology

- Frozen TX4 manuscript commit: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6` (`TX4_FINAL_MANUSCRIPT_FREEZE`).
- Latest independent validation implementation: committed at `069fa21204aa7c6f45fb44eb17e304ab99685a2b` on the prior validation branch before this clean worktree was created.
- Prior independent package retained separately under `research/ias2026_last_validation/`; its old GFL10 harness is not overwritten.
- Exact Python model source retained under `research/ias2026_last_validation/raw/true_same_model/canonical_source/models/`.
- Exact frozen GFL11 equations retained under `research/ias2026_last_validation/code/julia/frozen_gfl11.jl` and the prior source snapshots.
- Frozen flagship portfolio: `H4 = {30,33,35,37}`.
- Frozen P4 policy: `g=0.03625`, `k=1.425`, `t=1.5`, `h=1.0`.

## New branch

- Branch: `research/tx4-exact-p4-julia-reproduction`
- Parent: `069fa21204aa7c6f45fb44eb17e304ab99685a2b`
- Push: explicitly prohibited; no push will be performed.

## Expected gates

1. Exact GFL11 has 11 states per replaced converter.
2. H4 has `nx=86` (`4*11 + 6*7`).
3. All 16 P4 verdicts agree across Python and Julia.
4. H4 alpha difference is at most `1e-4 s^-1`.
5. H4 critical-frequency difference is at most `1e-4 Hz`.
6. Matched voltage-mode MAC is at least `0.95`.
7. Equilibrium voltage and angle differences are at most `1e-6`.
8. No hidden equation, policy, ordering, or state-count mismatch remains.

## Status at branch creation

No exact P4/GFL11 case has been executed on this branch before the preregistration commit.
