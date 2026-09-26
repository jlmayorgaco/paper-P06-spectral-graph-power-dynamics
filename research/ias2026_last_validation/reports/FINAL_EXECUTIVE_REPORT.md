> HISTORICAL PRE-LAST-VALIDATION SNAPSHOT. For the executed final verdict, use `LAST_INDEPENDENT_VALIDATION_EXECUTIVE_REPORT.md` and `../docs/CLAIM_LEDGER.md`.

# IAS 2026 Vancouver final scientific closure report

Date: 2026-09-19

## Executive decision

The final closure campaign does not clear the full award-level claim. It does
produce a defensible, reproducible audit package with explicit negative and
stopped results.

The strongest supported result is narrow: the frozen GFL11 device/port
implementation is numerically consistent across Python and Julia over 64
canonical cases and 602 transfer points, and its network harness satisfies the
declared residual limit. The official PowerDynamics and SimpleGFLDC IEEE-39
campaigns are valid alternative-model negative results. The official-model TDS
cases solve and provide decay diagnostics without a reliable modal-frequency
crossing.

The decisive limitation is that the exact Julia port of the frozen custom
SynchronousMachine/first-order AVR/two-state IEEEST PSS/no-governor model was
not completed. The true same-model cross-code gate is therefore
`STOPPED_BY_GATE`. Mechanism, physical robustness, genuine mixed holdout, and
new scaling claims are not promoted downstream.

## What survived

- G0 frozen-artifact reproduction: `FROZEN_ARTIFACT_REPRODUCED`.
- G1 canonical GFL11 Python/Julia device and transfer parity:
  `NUMERICALLY_VERIFIED`.
- G6 official PowerDynamics TDS execution: `NUMERICALLY_VERIFIED`, with scope
  limited to the official alternative model and decay diagnostics.
- G7 official SimpleGFLDC second-model census: `SECOND_MODEL_VALIDATED`, with
  a 16/16 stable negative result.
- G8 alternative-model holdout: `BLIND_HOLDOUT`, but single-class and not
  discriminative.
- The claim ledger and reviewer attack matrix provide safe wording and
  adversarial scope controls.

## What failed or remains unknown

- Exact same-model SG/AVR/PSS Python-to-Julia parity: `STOPPED_BY_GATE`.
- Full frozen IEEE-39 same-model base/V4 census: `STOPPED_BY_GATE`.
- Collective mechanism test: `NOT_TESTED` because its prerequisite stopped.
- Common physical uncertainty envelope, robust radius, and robust census:
  `NOT_TESTED`.
- New independent mixed holdout: not established; the executed holdout was
  single-class and needs a second preregistered holdout.
- New asymptotic scaling sweep: `NOT_TESTED`.

## Claim-safe conclusion

The package supports a scoped cross-code GFL11 parity result and documents
model-conditioned negative results. It does not support the stronger statement
that the frozen custom IEEE-39 collective blocker has been independently
validated, refuted, robustly generalized, or predictively detected.

## Artifact rule

This root contains no final paper or poster rewrite. The final ZIP contains
only this executive report as a PDF plus raw evidence, code, reports, docs,
figures, environment files, and hashes. Legacy generated paper/poster PDFs are
explicitly excluded.
