# P2 — Algebraic conditional retention frontier

Status: **PASS**, classified only as `EXACT_CONDITIONAL_RETENTION_ROOT`.

- Equation checked: `det(I₂+εN₃₈)=1+tr(N₃₈)ε+det(N₃₈)ε²` at `s=−0.05`, using the ExpN port current convention `Yinj=−Yraw`; roots were obtained from `−1/eig(N₃₈)`, then independently evaluated against the complete physical finite spectrum. A second root uses the historical numerical guard `δ_num=1.0e-9` s⁻¹.
- ExpN retained MW: `1.124344477653483`; guarded algebraic root: `1.1243444848621267` MW; absolute discrepancy `7.208643681977378e-9` MW. ε relative discrepancy: `6.411419025928842e-9`.
- Root α: exact `-0.049999997897093904` s⁻¹; guarded `-0.050000000241936744` s⁻¹. All physical poles pass: `true`.
- Ten single-SG architectures were tested algebraically; roots outside `[0,1]`, complex roots and roots dominated by another pole are retained/classified in the CSV.
- Scope: one-dimensional retention boundary at frozen ExpN gains. This is not a continuous free-gain optimum, robust optimum, or global certificate.
- New-design numerical guard audit: from the maximum observed ExpN analytic/PD alpha discrepancy `9.66594047857594e-10`, the proposed rule gives `δ_num=max(1e-6,10×error)=1e-6 s⁻¹`. The rank-two root has `ε_38=0.0013546412951520276`, retained `1.1243522749761832 MW` (increase `7.797322700175968e-6 MW`), and complete-spectrum alpha `−0.05000099857434908 s⁻¹`; all 101 finite physical poles were checked. This is a diagnostic root only: no KKT correction, new freeze, or PD validation was done, and the frozen candidate was not changed. See `TABLE_P2_guard_policy.csv` and `experiments/bnd_expP/audit_expP_numerical_guard.jl`.
- Time: `36.895` s. The algebraic root and full-pole gates permit P3. P3–P5 remain outstanding.
