# P2 — Algebraic conditional retention frontier

Status: **FAIL_ALGEBRAIC_ROOT**, classified only as `EXACT_CONDITIONAL_RETENTION_ROOT`.

- Equation checked: `det(I₂+εN₃₈)=1+tr(N₃₈)ε+det(N₃₈)ε²` at `s=−0.05`; roots were obtained from `−1/eig(N₃₈)`, then independently evaluated against the complete physical finite spectrum. A second root uses the historical numerical guard `δ_num=1.0e-9` s⁻¹.
- ExpN retained MW: `1.124344477653483`; guarded algebraic root: `23.73747364975028` MW; absolute discrepancy `22.613129172096798` MW. ε relative discrepancy: `20.112278417812483`.
- Root α: exact `-0.07896207490877126` s⁻¹; guarded `-0.0789620745994436` s⁻¹. All physical poles pass: `false`.
- Ten single-SG architectures were tested algebraically; roots outside `[0,1]`, complex roots and roots dominated by another pole are retained/classified in the CSV.
- Scope: one-dimensional retention boundary at frozen ExpN gains. This is not a continuous free-gain optimum, robust optimum, or global certificate.
- Time: `4.485` s. The algebraic root and full-pole gates permit P3. P3–P5 remain outstanding.
