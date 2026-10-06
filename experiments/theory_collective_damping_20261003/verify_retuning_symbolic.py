"""Symbolic verification of the theory extension; no grid experiments."""
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parent
checks = []

def check(name, expr):
    ok = sp.simplify(sp.expand_complex(expr)) == 0
    checks.append({"name": name, "passed": bool(ok),
                   "status": "EXACT_SYMBOLIC_ALGEBRA"})
    assert ok, name

alpha, h, tau = sp.symbols("alpha h tau", real=True)
omega, kp, ki = sp.symbols("omega kp ki", positive=True)
lam = alpha+sp.I*omega
conj_lam = alpha-sp.I*omega
kph = sp.exp(alpha*h)*(kp*sp.cos(omega*h)
                        +(ki+alpha*kp)*sp.sin(omega*h)/omega)
kih = sp.exp(alpha*h)*(ki*sp.cos(omega*h)
                        -((alpha**2+omega**2)*kp+alpha*ki)
                        *sp.sin(omega*h)/omega)
rotation = sp.exp(alpha*h)*(sp.cos(omega*h)+sp.I*sp.sin(omega*h))
check("finite_complex_PI_action_transport",
      lam*kph+kih-(lam*kp+ki)*rotation)
check("transport_recovers_baseline_Kp", kph.subs(h, 0)-kp)
check("transport_recovers_baseline_Ki", kih.subs(h, 0)-ki)
check("gain_flow_Kp", sp.diff(kph, h)-(2*alpha*kph+kih))
check("gain_flow_Ki", sp.diff(kih, h)+(alpha**2+omega**2)*kph)

Wr, Wi = sp.symbols("Wr Wi", real=True)
check("pole_assignment_gain_inverse",
      lam*Wi/omega+(Wr-alpha*Wi/omega)-(Wr+sp.I*Wi))

sr, si, tf, gr, gi = sp.symbols("sr si tf gr gi", real=True)
s = sr+sp.I*si
G = gr+sp.I*gi
F = s**2*(1+tf*s)-(s*kp+ki)*sp.exp(-s*tau)*G
check("gain_derivative_ratio_at_characteristic_operator",
      sp.diff(F, kp)-s*sp.diff(F, ki))
check("delay_gain_identity_at_characteristic_operator",
      sp.diff(F, tau)+s*(kp*sp.diff(F, kp)+ki*sp.diff(F, ki)))
# Strip the common nonzero exp(-s*tau) before checking the polynomial.
kernel_derivative = s*sp.diff(kph, h).subs(h, 0)+sp.diff(kih, h).subs(h, 0)-s*(s*kp+ki)
check("unavoidable_spillover_polynomial",
      kernel_derivative+kp*(s-lam)*(s-conj_lam))

mr, mi, cr, ci = sp.symbols("mr mi cr ci", real=True)
mu, gamma = mr+sp.I*mi, cr+sp.I*ci
mu_chain = -mu*(kp*mu*gamma+ki*gamma)+(ki+2*alpha*kp)*mu*gamma-(alpha**2+omega**2)*kp*gamma
check("full_modal_spillover_identity",
      mu_chain+kp*(mu-lam)*(mu-conj_lam)*gamma)
J = sp.Matrix([[sp.re(lam*gamma), sp.re(gamma)],
               [sp.im(lam*gamma), sp.im(gamma)]])
check("single_mode_gain_Jacobian_nonsingularity",
      J.det()+omega*(cr**2+ci**2))
check("Hopf_transport_Kp",
      kph.subs(alpha, 0)-(kp*sp.cos(omega*h)+ki*sp.sin(omega*h)/omega))
check("Hopf_transport_Ki",
      kih.subs(alpha, 0)-(ki*sp.cos(omega*h)-omega*kp*sp.sin(omega*h)))

result = {
    "all_passed": all(c["passed"] for c in checks),
    "check_count": len(checks), "checks": checks,
    "simulation_count": 0, "optimization_runs": 0,
    "spectral_sweeps": 0,
    "interpretation": "Algebra verification only; no IEEE-39 claim or novelty certificate",
}
(OUT/"RETUNING_SYMBOLIC_CHECKS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
prior = json.loads((OUT/"PROVENANCE.json").read_text(encoding="utf-8"))
provenance = {
    "revision": "theory-only retuning extension, 2026-10-03",
    "python": platform.python_version(), "sympy": sp.__version__,
    "parent_provenance_sha256": sha(OUT/"PROVENANCE.json"),
    "parent_theory_sha256": next(e["sha256"] for e in prior["new_artifacts"]
                                if e["path"] == "THEORY.tex"),
    "current_theory_sha256": sha(OUT/"THEORY.tex"),
    "verification_script_sha256": sha(Path(__file__)),
    "new_checks_sha256": sha(OUT/"RETUNING_SYMBOLIC_CHECKS.json"),
    "prior_results_preserved": True,
    "scope": "No trajectory, root search, optimization, commit, or push",
}
(OUT/"RETUNING_PROVENANCE.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
print(json.dumps({k: result[k] for k in ["all_passed", "check_count", "simulation_count",
                                        "optimization_runs", "spectral_sweeps"]}, indent=2))
