"""Exact algebra only: no trajectories, optimization, or spectral sweeps."""
import hashlib
import json
import platform
import subprocess
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
checks = []

def check(name, expression):
    if isinstance(expression, sp.MatrixBase):
        ok = all(sp.simplify(e) == 0 for e in expression)
    else:
        ok = sp.simplify(expression) == 0
    checks.append({"name": name, "passed": bool(ok),
                   "scope": "EXACT_SYMBOLIC_ALGEBRA_NOT_DYNAMIC_VALIDATION"})
    assert ok, name

s = sp.symbols("s")
V, kp, ki, tau, tf = sp.symbols("V kp ki tau tf", positive=True)
H = V*(kp*s+ki)*sp.exp(-s*tau)/(s**2*(1+tf*s)+V*(kp*s+ki)*sp.exp(-s*tau))
expected = 1-s**2/(V*ki)+(kp/ki-tf-tau)*s**3/(V*ki)
check("actual_PLL_series_through_cubic", sp.series(H, s, 0, 4).removeO()-expected)

# Finite Schur update; rational example with symbolic update size.
kappa = sp.symbols("kappa")
D = sp.Matrix([[2, 1, 0, 1], [0, 3, 1, -1],
               [1, 0, 4, 1], [-1, 2, 0, 5]])
b = sp.Matrix([1, 2, -1, 3])
c = sp.Matrix([2, -1, 1, 1])
M = sp.diag(2, 3)
def schur(X):
    return M*(X[:2, :2]-X[:2, 2:]*X[2:, 2:].inv()*X[2:, :2])
hidden_inv = D[2:, 2:].inv()
u = M*(b[:2, :]-D[:2, 2:]*hidden_inv*b[2:, :])
vrow = c[:2, :].T-c[2:, :].T*hidden_inv*D[2:, :2]
h = (c[2:, :].T*hidden_inv*b[2:, :])[0]
check("finite_physical_Schur_update",
      schur(D+kappa*b*c.T)-schur(D)-kappa/(1+kappa*h)*u*vrow)

# Canonical two-dimensional span used in the general Hermitian proof.
A, br, bi, C, z = sp.symbols("A br bi C z", real=True)
a = sp.Matrix([A, 0])
v = sp.Matrix([br+sp.I*bi, C])
Hrank = (a*v.conjugate().T+v*a.conjugate().T)/2
char = Hrank.charpoly(z)
# charpoly creates its own generator without the supplied real assumption.
char_expression = char.as_expr().subs(char.gen, z)
check("Hermitian_rank_two_characteristic",
      char_expression-(z**2-A*br*z-A**2*C**2/4))
check("Hermitian_eigenvalue_product",
      Hrank.det()+((a.conjugate().T*a)[0]*(v.conjugate().T*v)[0]
                   -sp.Abs((v.conjugate().T*a)[0])**2)/4)

# General two-node check includes heterogeneous inertia and nonsymmetric D0.
m1, m2, d11, d12, d21, d22, w = sp.symbols(
    "m1 m2 d11 d12 d21 d22 w", positive=True)
mass = sp.diag(m1, m2)
one = sp.ones(2, 1)
m = m1+m2
damp = sp.Matrix([[d11, d12], [d21, d22]])
L = w*sp.Matrix([[1, -1], [-1, 1]])
Ldag = sp.Matrix([[1, -1], [-1, 1]])/(4*w)
poly = s**2*mass+s*damp+L
row = one.T*mass/m
column = mass*one/m
transfer = sp.cancel(s*(row*poly.adjugate()*column)[0]/poly.det())
zcoi = sp.cancel(1/transfer)
d0 = (one.T*damp*one)[0]
bR = damp*one-d0*mass*one/m
bL = damp.T*one-d0*mass*one/m
mapp = m-(bL.T*Ldag*bR)[0]
check("general_two_node_COI_DC", zcoi.subs(s, 0)-d0)
check("general_two_node_COI_moment",
      sp.diff(zcoi, s).subs(s, 0)-mapp)
v1, v2 = sp.symbols("v1 v2", real=True)
vinput = sp.Matrix([v1, v2])
P = v1+v2
Hv = sp.cancel(s*(row*poly.adjugate()*vinput)[0]/poly.det())
area = -sp.diff(Hv, s).subs(s, 0)
ainput = vinput-P*mass*one/m
check("general_two_node_input_location_area",
      area-P*mapp/d0**2-(bL.T*Ldag*ainput)[0]/d0)

example = {m1: 1, m2: 1, d11: 1, d12: 0, d21: 0, d22: 9, w: 1}
Gexample = (s**2+5*s+2)/(2*(s**3+10*s**2+11*s+10))
check("stable_overshoot_example_transfer", transfer.subs(example)-Gexample)
check("stable_overshoot_example_moment", mapp.subs(example)+14)
check("stable_overshoot_example_signed_area",
      -sp.diff(Gexample, s).subs(s, 0)+sp.Rational(14, 100))
check("homogeneous_example_is_monotone_first_order",
      transfer.subs({**example, d11: 5, d22: 5})-1/(2*s+10))
assert 10*11 > 10  # Routh-Hurwitz, with all cubic coefficients positive.

result = {
    "method": "symbolic exact algebra; written proofs are in THEORY.tex",
    "simulation_count": 0, "optimizer_runs": 0,
    "checks": checks,
    "all_passed": all(c["passed"] for c in checks),
    "example": {"network": "illustrative two-node, not IEEE-39",
                "m": 2, "d0": 10, "graph_correction": 16, "m_app": -14,
                "signed_area_per_unit_step": "-7/50",
                "stable_by": "Routh-Hurwitz: 10*11 > 10"}}
(OUT/"SYMBOLIC_CHECKS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

paths = [
    "experiments/nonlinear_codesign_20261001/ReducedDAE.jl",
    "src/bnd_model_expN/PDExactDesignN.jl",
    "experiments/graph_gsp_codesign_20261003/DelayedEvents.jl",
    "experiments/physical_collective_damping_20261003/THEORY.md",
    "experiments/physical_collective_damping_20261003/action_attribution.py",
    "experiments/graph_gsp_codesign_20261003/THEORY.md",
]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
provenance = {
    "date": "2026-10-03",
    "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "git_status_porcelain": subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
    "python": platform.python_version(), "sympy": sp.__version__,
    "source_files": [{"path": p, "sha256": sha(ROOT/p)} for p in paths],
    "new_artifacts": [{"path": p.name, "sha256": sha(p)}
                      for p in [OUT/"THEORY.tex", Path(__file__), OUT/"SYMBOLIC_CHECKS.json"]],
    "scope": "theory and algebra only; no frozen file modified; no commit/push",
    "full_model_claim": "no new IEEE-39 feasibility/optimum/nonlinear validation",
}
(OUT/"PROVENANCE.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
print(json.dumps({"all_passed": result["all_passed"], "check_count": len(checks),
                  "simulation_count": 0, "output": str(OUT)}, indent=2))
