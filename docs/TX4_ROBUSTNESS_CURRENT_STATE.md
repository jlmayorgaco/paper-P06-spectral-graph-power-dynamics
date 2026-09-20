# TX4 Robustness Campaign — Current State Before New Results

Date: 2026-09-20

## Frozen parent and branch discipline

- Frozen exact P4/GFL11 closure parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`.
- New worktree branch: `research/tx4-final-robustness-statistics`.
- Frozen exact worktree is a separate sibling and is not modified by this campaign.
- No push is authorized.

## Archaeology result

The final-closure dashboard records G10 (common physical uncertainty envelope),
G11 (robust V4 census), and G12 (g-by-epsilon atlas) as `NOT_TESTED`. The
retained uncertainty manifest is a specification-only declaration with no
weights and explicitly forbids portfolio-specific refitting and calling a
closure-space radius a physical radius. The archived robustness report also
states that no physical robust radius, common envelope, robust V4 census, or
g-by-epsilon atlas was completed.

The exact parent does contain the frozen same-model IEEE-39 Python/Julia P4
implementation and the exact GFL11 H4 closure. It also contains the prior
modal-scope separation table, including the historical g-boundary witness
`g*=0.2076814051903784`, but that table is retrospective evidence and is not
treated as a new robustness result.

## Scope carried forward

The new campaign uses the frozen reduced semi-explicit IEEE-39 DAE and the
matched-reactive-power policy. It reports nominal-closure metrics separately
from robustness metrics and does not promote results to EMT, switching,
current-limit, DC-link, protection, hardware, universal IEEE-39, or
all-model claims.

## Known dependency limitation

The prior P1–P6 audit stopped full PowerDynamics same-model network parity at
the dynamic infinite-bus gate. Therefore the new campaign is a Python exact
reduced-DAE robustness campaign with a separately reported 32-condition Julia
cross-code spot-check; it is not a fresh PowerDynamics parity certificate.

## Pre-existing numerical probe

A single unrecorded timing probe was used only to size the execution plan
before this preregistration was committed. It is excluded from every result
file and claim. This deviation is recorded in
`docs/TX4_ROBUSTNESS_DEVIATIONS.md`.
