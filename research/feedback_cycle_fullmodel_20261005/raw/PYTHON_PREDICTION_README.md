# Python pre-Julia closure: feedback-cycle resonance and sparse repair

## Decision
**GO to PowerDynamics/Julia as a frozen falsification test.**

The Python campaign now provides a specific prediction, rather than an open-ended search.

## 1. Nonlinear 9-bus DDE: resonance is real in time domain
The existing 21-ms model has a rightmost controller mode at approximately
`+0.88615 ± j 48.6485 s^-1` (`7.74265 Hz`). At that root the locally dressed
interaction matrix has an eigenvalue at `-1` while the three local factors stay
nonsingular, so the root is a collective closure rather than a local singularity.

With a 0.5-MW sinusoidal load at the modal frequency:
- 20-ms stable baseline: max PLL detector error ≈ **0.505°**
- 21-ms no retuning: ≈ **32.435°**
- 21-ms sparse repair: ≈ **0.678°**
- 21-ms analytical transport: ≈ **0.485°**

The same unstable case has a much smaller SG-frequency excursion than its PLL-state
growth, so this is primarily a synchronization/control resonance rather than an
ordinary electromechanical swing.

## 2. Exploratory reduced IEEE-39
This is **not** PowerDynamics parity. It is a new Python mechanism model assembled
from the public solved JuliaEnergy/IEEE39.jl static case and its machine H/Xd' data,
with reduced SG/GFL dynamics. It exists only to make a falsifiable prediction.

Main results:
- all-SG rightmost oscillatory family: roughly **1–1.6 Hz**
- after 75% GFL sharing at buses 30–38: a distinct **7.5–7.9 Hz** PLL/controller family
- tracked branch crosses the imaginary axis at
  **20.3025018 ms, 7.8792356 Hz**
- at the crossing, the **35↔36** two-cycle is almost closed:
  `|q35,36 q36,35| = 0.9961304`, phase `−0.2582°`, `|1-p| = 0.0059335`
- other important physical/dynamic pairs are **30↔37** and **33↔34**
- graph/Kron screening plus modal gain authority selects buses **30,33,36,37**
- a small exact-root iterative corrector moves the rightmost root
  from approximately **+0.97063** to **−0.050008 s^-1**
- restoring all Kp to baseline and retaining only the Ki changes at those four sites
  still gives approximately **−0.051786 s^-1**
- the Kp-only version fails; this is consistent with the separate negative IEEE-39
  Kp-only experiment reported by Claude
- removing any one of the four Ki changes from this constructed repair loses the
  declared −0.05 margin; this is not a proof that no alternative three-site design exists

## 3. What is worth testing in Julia
Freeze these predictions before looking at the PowerDynamics result:
1. all-SG model should not contain the same 7–8 Hz PLL family;
2. SG→GFL replacement should introduce a controller/synchronization family;
3. the limiting full-model feedback core should involve the same electrical neighborhoods
   (35–36, 30–37, 33–34) or decisively falsify the reduced prediction;
4. Ki sensitivity should dominate Kp sensitivity for the critical family;
5. the preregistered candidate support `{30,33,36,37}` should be tested before running
   any unrestricted optimizer;
6. forcing near the limiting frequency should excite PLL/control states more strongly
   than matched off-resonance forcing;
7. full root counts and nonlinear event guards remain the final decision.

## Scope
The numerical IEEE-39 threshold `20.3025 ms` is model-specific and is not a physical
universal latency limit. The reduced Python IEEE-39 does not include full
PowerDynamics SimpleGFLDC/AVR/PSS fidelity. The 9-bus nonlinear result is stronger
time-domain evidence but is a separate declared reduced model.

The key result at this stage is therefore **mechanistic and predictive**:
the network can close frequency-specific feedback cores among local GFL controllers,
and graph/core screening can greatly reduce the sparse retuning search space.
