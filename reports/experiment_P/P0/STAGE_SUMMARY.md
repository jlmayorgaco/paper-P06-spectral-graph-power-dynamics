# P0 — ExpN reproduction and P/Q contract

Status: **PASS**. This is an independent rebuild of the frozen ExpN candidate and one deterministic withheld mixed point.

- Equation checked: frozen explicit SG/GFL P/Q split at the original bus voltage, with the original bus-31 and bus-39 ZIP loads retained component-wise; physical descriptor and quotient spectrum from ExpN.
- Maximum trim residual: `5.395296582857102e-12`; P/Q errors: `3.5583980206865816e-13` / `4.632738637155853e-14` pu.
- Reduced-A relative error: `2.644242322428162e-16`; complete physical pole matching error: `4.2772862799896e-8` s⁻¹.
- Candidate: retained `1.124344477653483` MW, converted `5401.636745501193` MW, α=`-0.0500000011050941` s⁻¹. Candidate hash remains `3915a8f57da0f552779b11464572f14a86f0b648b532be77124f668eeff9f263`.
- Certificate scope: **independent ExpN regression and two point identity check only**; no new optimum claim.
- Time: `285.52` s; analytic design PowerDynamics calls: 0; independent validation builds: 2.
- This permits P1 only because both P/Q/spectrum identity gates pass. No prior artifact was modified.
