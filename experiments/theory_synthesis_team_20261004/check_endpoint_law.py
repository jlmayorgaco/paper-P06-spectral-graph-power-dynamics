"""Algebra and static-endpoint checks only. No dynamical simulation.

The example is dimensionless and is NOT the IEEE-39 PLL model.
Run from anywhere: python <this file>
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
P = sp.Integer(2)
m = sp.Rational(41, 40)
d0 = sp.Integer(4)
delta = -sp.pi / 6
E_sec = sp.pi / 3
area = sp.simplify(P * (m - E_sec) / d0**2)
area_linear = sp.simplify(P * (m - 1) / d0**2)
energy_floor = sp.simplify(-d0 * area)
assert sp.simplify(area - (123 - 40 * sp.pi) / 960) == 0
assert area_linear == sp.Rational(1, 320)
assert float(area) < 0 < float(area_linear)
assert sp.simplify(sp.sin(delta) + sp.Rational(1, 2)) == 0
assert float(1 - sp.pi / 4 - m / 8) > 0

# Independent symbolic balance: mass/damping decomposition implies the area law.
mass, damp, omega, bdelta, Q, kappa = sp.symbols("m d Omega bdelta Q kappa", nonzero=True)
A = sp.symbols("A")
balance = (mass + kappa) * omega - damp * A + bdelta - Q
solved_A = sp.solve(balance, A)[0]
assert sp.simplify(solved_A - ((mass + kappa) * omega + bdelta - Q) / damp) == 0

# Nonlinear six-node endpoint pair. No small-angle substitution.
edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 3), (1, 4)]
B = np.zeros((6, len(edges)))
for e, (i, j) in enumerate(edges):
    B[i, e], B[j, e] = 1, -1
gamma = np.array([2., 3., 1.5, 2.5, 4., 2.2, 1.7, 2.8])
q0 = np.array([.08, -.12, .09, -.07, .03, -.01])
q1 = np.array([-.04, .13, -.08, .11, -.02, .04])
eta0, eta1 = B.T @ q0, B.T @ q1
delta_q = q1 - q0
DeltaF = B @ (gamma * (np.sin(eta1) - np.sin(eta0)))
Om = .25
b = -DeltaF / Om
weights = gamma * (np.sin(eta1) - np.sin(eta0)) / (eta1 - eta0)
Lbar = (B * weights) @ B.T
pinv = np.linalg.pinv(Lbar, hermitian=True)
potential = pinv @ b
E = float(b @ potential)
flow = weights * (B.T @ potential)
flow_energy = float(np.sum(flow**2 / weights))
E_incremental = float(delta_q @ DeltaF)
secant_error = float(np.linalg.norm(Lbar @ delta_q - DeltaF))
assert secant_error < 1e-12
assert abs(E - flow_energy) < 1e-11
assert abs(E - E_incremental / Om**2) < 1e-11
vals, vecs = np.linalg.eigh(Lbar)
spectral_energy = float(np.sum((vecs[:, 1:].T @ b)**2 / vals[1:]))
assert abs(E - spectral_energy) < 1e-11
cut_bounds = []
for mask in range(1, 2**6 - 1):
    ix = np.array([(mask >> j) & 1 for j in range(6)], dtype=bool)
    crossing = np.array([ix[i] != ix[j] for i, j in edges])
    cut_bounds.append(float(b[ix].sum()**2 / weights[crossing].sum()))
assert max(cut_bounds) <= E + 1e-11
mass_vec = np.array([1.2, .8, 1.1, .7, 1.3, .9])
dc_damping_total = 2 * max(1., float(np.max(-b * mass_vec.sum() / mass_vec)))
damping_vec = dc_damping_total * mass_vec / mass_vec.sum() + b
assert np.min(damping_vec) > 0
assert np.linalg.norm(DeltaF - (Om * dc_damping_total * mass_vec / mass_vec.sum() - Om * damping_vec)) < 1e-12

result = {
    "status": "PROVED_REDUCED_MODEL / NUMERICALLY_VALIDATED_ALGEBRA_ONLY",
    "dynamic_simulations_executed": 0,
    "model": "Dimensionless nonlinear lossless swing; not IEEE-39 or a PLL realization",
    "example": {
        "M_diagonal": [str(sp.Rational(41, 80))] * 2,
        "D_diagonal": [3, 1], "line_sine_coefficient": 1,
        "step": [1, 1], "final_frequency": .5,
        "final_angle_difference": str(delta),
        "secant_mismatch": str(E_sec),
        "signed_frequency_area_exact": str(area),
        "signed_frequency_area_numeric": float(area),
        "signed_frequency_area_linearized": str(area_linear),
        "net_opposing_energy_strict_lower_bound": str(energy_floor),
        "net_opposing_energy_lower_bound_numeric": float(energy_floor),
        "initial_rotating_frame_energy": str(m / 8),
        "lowest_cohesive_boundary_energy": str(1 - sp.pi / 4),
        "energy_barrier_margin": float(1 - sp.pi / 4 - m / 8),
        "claim": "Stable convergent baseline with unavoidable endpoint overshoot; the linear signed-area test misses this case."
    },
    "six_node_static_check": {
        "secant_identity_residual": secant_error,
        "E_pseudoinverse": E, "E_flow": flow_energy, "E_gsp": spectral_energy,
        "E_incremental_divided_by_Omega_squared": E_incremental / Om**2,
        "max_cut_lower_bound": max(cut_bounds),
        "positive_damping_min": float(damping_vec.min()),
        "note": "Endpoint algebra only; no assertion of dynamical convergence for this constructed pair."
    },
    "versions": {"python": platform.python_version(), "numpy": np.__version__, "sympy": sp.__version__},
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(ROOT / "ENDPOINT_LAW_CHECKS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

# Static analytical curve for the same two-node family, fixed total mass.
p = np.linspace(.001, 3.9, 400)
E_curve = 4 * np.arcsin(p / 4) / p
area_curve = p / 16 * (float(m) - E_curve)
linear_curve = p / 16 * (float(m) - 1)
fig, ax = plt.subplots(figsize=(8.1, 4.7), constrained_layout=True)
ax.plot(p, area_curve, label="Exact nonlinear endpoint law", color="#116466", lw=2.5)
ax.plot(p, linear_curve, label="Initial-Jacobian moment", color="#b36d16", lw=2, ls="--")
ax.axhline(0, color=".4", lw=.8)
ax.scatter([2], [float(area)], color="#116466", zorder=4)
ax.annotate("Proved convergent example\nA < 0: overshoot is unavoidable", (2, float(area)),
            xytext=(.18, -.041), arrowprops={"arrowstyle": "->", "color": ".3"}, fontsize=10)
ax.set(xlabel="Step amplitude P (dimensionless)", ylabel="Signed frequency area A (dimensionless)",
       title="Two equilibria reveal a limitation missed by the initial Jacobian")
ax.text(.02, .04, "Analytical endpoints only. Convergence proved at P = 2.\nNot IEEE-39 data; not a replacement-capacity frontier.",
        transform=ax.transAxes, fontsize=9, color=".3")
ax.legend(loc="lower left", bbox_to_anchor=(0, .17), frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(alpha=.18)
fig.savefig(ROOT / "FIG_ENDPOINT_AREA_ANALYTIC.png", dpi=220)
fig.savefig(ROOT / "FIG_ENDPOINT_AREA_ANALYTIC.svg")
plt.close(fig)
print(json.dumps({"area": float(area), "energy_lower_bound": float(energy_floor),
                  "algebra_residual": secant_error, "dynamic_simulations": 0}, indent=2))
