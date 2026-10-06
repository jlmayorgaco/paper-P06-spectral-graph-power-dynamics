"""Exact algebra checks for the design closure. No power-grid experiments."""
import hashlib
import json
import platform
from pathlib import Path

import sympy as sp

OUT = Path(__file__).resolve().parent
checks = []
PRESERVED = {
    "PROVENANCE.json": "cfa4e6c752688afb86c1aee364e9dd46405d1d2677ac6bbd19cc331f923dd50e",
    "SYMBOLIC_CHECKS.json": "f9a50900de348efaf0a72f0927d64fb7e5e8eefadff404893c1efeacfb9e00fe",
    "RETUNING_PROVENANCE.json": "6848146ee623885f54b43e3cea09990563019b3e8f5fbbae77aa57a5655ce4b0",
    "RETUNING_SYMBOLIC_CHECKS.json": "9c355fa00655daade70efa3d759ee00f284e93dbcc5d604a625e9820217dc59b",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(name, ok):
    checks.append({"name": name, "passed": bool(ok),
                   "status": "EXACT_IDENTITY", "verification": "symbolic algebra"})
    assert ok, name


def zero(name, expression):
    values = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    record(name, all(sp.simplify(v) == 0 for v in values))


alpha, omega, wr, wi = sp.symbols("alpha omega wr wi", real=True, nonzero=True)
lam = alpha + sp.I * omega
kp, ki = wi / omega, wr - alpha * wi / omega
zero("unique_real_gains", lam * kp + ki - (wr + sp.I * wi))
C1 = sp.Matrix([[alpha, 1], [omega, 0]])
zero("single_pair_inverse", C1 * sp.Matrix([kp, ki]) - sp.Matrix([wr, wi]))
record("single_pair_rank", sp.simplify(C1.det() + omega) == 0)

# Exact two-site row reconstruction at arbitrary complex characteristic samples.
# These samples are algebraic fixtures, not a physical grid or stability result.
l = -1 + 2 * sp.I
q = sp.Matrix([1, 2 + sp.I])
p0 = sp.Matrix([2, 3])
i0 = sp.Matrix([4, 5])
tf = [sp.Rational(1, 10), sp.Rational(1, 20)]
d1, d2 = sp.symbols("d1 d2", nonzero=True)
delay = sp.diag(d1, d2)
D = sp.diag(*[l ** 2 * (1 + f * l) for f in tf])
y = sp.Matrix([D[j, j] * q[j] / (delay[j, j] * (l * p0[j] + i0[j]))
               for j in range(2)])
zero("full_matrix_row_assignment",
     D * q - (l * sp.diag(*p0) + sp.diag(*i0)) * delay * y)
for j in range(2):
    W = sp.simplify(D[j, j] * q[j] / (delay[j, j] * y[j]))
    zero(f"site_{j}_gain_recovery_p", sp.im(W) / sp.im(l) - p0[j])
    zero(f"site_{j}_gain_recovery_i", sp.re(W) - sp.re(l) * p0[j] - i0[j])
c, qv, yv, sv, tv, fv = sp.symbols("c q y s tau tf", nonzero=True)
W = sv ** 2 * (1 + fv * sv) * sp.exp(sv * tv) * qv / yv
zero("complex_pattern_scale_invariance", W.subs({qv: c*qv, yv: c*yv}) - W)

# Total derivative with variable pole, pattern and plant response; tau fixed.
t = sp.symbols("t", real=True)
st, qt, yt = (sp.Function(n)(t) for n in ("s", "q", "y"))
Wt = st**2 * (1 + fv*st) * sp.exp(st*tv) * qt / yt
dWt = Wt*((2/st + fv/(1+fv*st) + tv)*sp.diff(st, t)
          + sp.diff(qt, t)/qt - sp.diff(yt, t)/yt)
zero("total_gain_quotient_derivative", sp.diff(Wt, t) - dWt)
ar, om, reW, imW = (sp.Function(n)(t) for n in ("alpha", "omega", "Wr", "Wi"))
kpt = imW / om
kit = reW - ar * kpt
dkpt = (sp.diff(imW, t) - kpt * sp.diff(om, t)) / om
zero("real_proportional_gain_chain_rule", sp.diff(kpt, t) - dkpt)
zero("real_integral_gain_chain_rule",
     sp.diff(kit, t) - (sp.diff(reW, t)-ar*dkpt-kpt*sp.diff(ar, t)))

# Two modes, one local signal and one proposed neighbor. Exact rational fixture.
# lambda = -1+i, -2+3i; local samples 1,1; neighbor samples i,2+i.
poles = [-1 + sp.I, -2 + 3*sp.I]
def columns(samples):
    rows = []
    for s, z in zip(poles, samples):
        rows.extend([[sp.re(s*z), sp.re(z)], [sp.im(s*z), sp.im(z)]])
    return sp.Matrix(rows)
A = columns([1, 1])
B = columns([sp.I, 2+sp.I])
x0 = sp.Matrix([2, 5])
u0 = sp.Matrix([1, 0])
b = A*x0 + B*u0
Pi = sp.eye(4) - A*A.pinv()
r = Pi*b
H = Pi*B
record("local_two_pair_conflict_is_nonzero", r != sp.zeros(4, 1))
record("neighbor_adds_the_missing_rank", A.row_join(B).rank() == 4)
zero("projector_annihilates_self_gain_actions", Pi*A)
zero("projector_is_idempotent", Pi*Pi-Pi)
R = sp.diag(2, 3)
Gram = H*R.inv()*H.T
ustar = R.inv()*H.T*Gram.pinv()*r
zero("added_gain_formula_satisfies_conflict", H*ustar-r)
zero("added_gain_formula_recovers_unique_solution", ustar-u0)
xstar = A.pinv()*(b-B*ustar)
zero("self_gain_reconstruction", xstar-x0)
zero("full_assignment_after_neighbor_correction", A*xstar+B*ustar-b)
zero("minimum_added_gain_cost", (ustar.T*R*ustar)[0]-(r.T*Gram.pinv()*r)[0])
zero("weighted_stationarity", R*ustar-H.T*Gram.pinv()*r)

# Compatible targets do not create an artificial conflict.
zero("compatible_local_demands", Pi*(A*x0))
# Rank failure can persist despite adding a neighbor with redundant signals.
B_redundant = 2*A
H_redundant = Pi*B_redundant
zero("redundant_neighbor_has_no_added_authority", H_redundant)
record("redundant_neighbor_cannot_resolve_conflict",
       H_redundant.row_join(r).rank() > H_redundant.rank())

w, a0, a1, a2, e1, e2 = sp.symbols("w a0 a1 a2 e1 e2", real=True)
L = sp.Matrix([[w, -w], [-w, w]])
K = a0*sp.eye(2) + a1*L + a2*(L**2)
E = sp.diag(e1, e2)
zero("polynomial_delay_commutator", L*K*E-K*E*L-K*(L*E-E*L))
phase1, phase2 = sp.symbols("phase1 phase2", real=True)
entry_norm_sq = ((sp.cos(phase1)-sp.cos(phase2))**2
                 +(sp.sin(phase1)-sp.sin(phase2))**2)
zero("exact_frequency_edge_identity",
     sp.trigsimp(2*w**2*entry_norm_sq-8*w**2*sp.sin((phase1-phase2)/2)**2))

for name, expected in PRESERVED.items():
    assert sha(OUT/name) == expected, f"Prior record changed: {name}"
root = OUT.parents[1]
prior_sources = json.loads((OUT/"PROVENANCE.json").read_text(encoding="utf-8"))["source_files"]
for source in prior_sources:
    assert sha(root/source["path"]) == source["sha256"], f"Prior source changed: {source['path']}"

result = {
    "all_passed": all(c["passed"] for c in checks),
    "check_count": len(checks), "checks": checks,
    "simulation_count": 0, "optimization_runs": 0, "spectral_sweeps": 0,
    "interpretation": "Proof supplements and algebraic fixtures only; no new IEEE-39 design",
}
(OUT/"CLOSURE_SYMBOLIC_CHECKS.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
provenance = {
    "revision": "analytical modal compatibility and design closure, 2026-10-03",
    "parent_theory_sha256": "9d7469fecb1692eb5063a5dc8e4624128b062b979fc99ba413a366175a47da5d",
    "current_theory_sha256": sha(OUT/"THEORY.tex"),
    "verification_script_sha256": sha(Path(__file__)),
    "checks_sha256": sha(OUT/"CLOSURE_SYMBOLIC_CHECKS.json"),
    "python": platform.python_version(), "sympy": sp.__version__,
    "preserved_prior_record_hashes": PRESERVED,
    "preserved_model_source_hashes": prior_sources,
    "scope": "No grid simulation, optimization, spectral run, commit or push",
}
(OUT/"CLOSURE_PROVENANCE.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8")
print(json.dumps({"all_passed": result["all_passed"], "check_count": len(checks),
                  "grid_simulations": 0, "optimizations": 0}))
