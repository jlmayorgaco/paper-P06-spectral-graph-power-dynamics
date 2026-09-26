# P1 failure-mode record

This file preserves the integration failures encountered while constructing
the one-device/infinite-bus PowerDynamics harness. The independent device
equation check passed in every attempt.

| attempt | change | result |
|---|---|---|
| 1 | Campaign Project.toml omitted the direct `SciCompDSL` dependency. | `Package SciCompDSL not found`; fixed by adding the existing manifest-pinned UUID to the project. |
| 2 | Three compiled bus objects had no unique graph names. | `Either all vertex models must have assigned graphelement or all vertex models must have unique name`; fixed with `@named`. |
| 3 | Current-source GFL plus junction plus micro-line plus slack. | Dynamic initialization residual was `1.96e-7` at the requested default tolerance; relaxing the harness tolerance exposed the next topology issue. |
| 4 | Current-source GFL loopback directly to a slack bus. | The returned network state had zero retained states; the residual reduction therefore could not validate the dynamic device. |

The result is a clean P1 integration stop, not a numerical parity mismatch.
The exact equation oracle, Julia transcription, compile result, equilibrium
residual, and transfer sweep remain in the adjacent raw files.
