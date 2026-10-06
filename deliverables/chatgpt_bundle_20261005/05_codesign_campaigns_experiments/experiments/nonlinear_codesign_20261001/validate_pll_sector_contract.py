"""Independent nonlinear consistency/falsification checks of a verified PLL contract.

The rational matrix proof supplies the all-input guarantee. These finite tests
only check equation/sign/coordinate and implementation consistency.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/pll_sector_contract"
data = json.loads((OUT/"synthesis.json").read_text(encoding="utf-8"))
cert = json.loads((OUT/"verification.json").read_text(encoding="utf-8"))
assert cert["synthesis_sha256"] == hashlib.sha256((OUT/"synthesis.json").read_bytes()).hexdigest()
p, best = data["specification"], data["best"]
kp, ki = best["Kp"], best["Ki"]
eps = cert["certificates"][-1]["epsilon_verified"]
X = np.array(best["X_scaled"])
P = np.linalg.inv(X)
S = np.diag([p["angle_max"]]+[2*np.pi*p["frequency_max_hz"]]*3)
Sinv = np.linalg.inv(S)
L = np.linalg.cholesky(X)
rng = np.random.default_rng(20261001)
tau, taum, alpha, kappa = [p[k] for k in ("pll_tau", "measurement_tau", "decay_rate", "frame_rate_per_voltage")]


def vectorfield(z, V0, w):
    delta, omega, xi, chi = z
    ud, uq, nu = eps*w[0], eps*w[1], kappa*eps*w[2]
    eu = -np.sin(delta)*ud+np.cos(delta)*uq
    error = -V0*np.sin(delta)+eu
    f = np.array([omega-nu, (xi+kp*error-omega)/tau,
                  ki*error, (omega-chi)/taum])
    return f, np.array([eu/eps, w[2]])


max_residual = -np.inf
max_boundary_derivative = -np.inf
for case in range(5000):
    direction = rng.normal(size=4); direction /= np.linalg.norm(direction)
    radius = 1.0 if case % 2 else rng.random()**0.25
    y = radius*(L@direction)
    w = rng.normal(size=3); w /= np.linalg.norm(w)
    V0 = rng.uniform(p["voltage_min"], p["voltage_max"])
    f, projected_input = vectorfield(S@y, V0, w)
    vv, vd = y@P@y, 2*y@P@Sinv@f
    max_residual = max(max_residual, vd+2*alpha*vv-2*alpha*(projected_input@projected_input))
    if radius == 1.0:
        max_boundary_derivative = max(max_boundary_derivative, vd)
assert max_residual < 1e-7 and max_boundary_derivative < 1e-7


def adversary(y):
    delta = (S@y)[0]
    bb = np.array([0, kp/tau, ki, 0])
    bn = np.array([-1, 0, 0, 0])
    forcing = 2*y@P@Sinv@bb*eps
    direction = np.array([-forcing*np.sin(delta), forcing*np.cos(delta),
                          2*y@P@Sinv@bn*kappa*eps])
    norm = np.linalg.norm(direction)
    return direction/norm if norm > 1e-14 else np.array([0.0, 1.0, 0.0])


cR = np.array([0, p["frequency_max_hz"]/(taum*p["rocof_max_hz_s"]),
               0, -p["frequency_max_hz"]/(taum*p["rocof_max_hz_s"])])
y_boundary = X@cR/np.sqrt(cR@X@cR)
trajectories = []
histories = {}
for V0 in (p["voltage_min"], p["voltage_max"]):
    for initial_name, initial in (("origin", np.zeros(4)), ("rocof_boundary", y_boundary)):
        for input_name in ("chirp", "state_adversary"):
            def rhs(t, y):
                if input_name == "state_adversary":
                    w = adversary(y)
                else:
                    w = np.array([np.sin(1.0+t*t), np.sin(7*t+0.4), np.cos(3*t)])
                    w /= max(1.0, np.linalg.norm(w))
                return Sinv@vectorfield(S@y, V0, w)[0]
            sol = solve_ivp(rhs, (0, 3.0), initial, method="Radau", rtol=2e-8,
                            atol=1e-10, dense_output=True, max_step=0.01)
            assert sol.success
            t = np.linspace(0, 3, 3001)
            yy = sol.sol(t)
            zz = S@yy
            storage = np.einsum("ij,ji->i", yy.T@P, yy)
            freq = zz[1]/(2*np.pi)
            rocof = (zz[1]-zz[3])/(2*np.pi*taum)
            key = f"V{V0}_{initial_name}_{input_name}"
            metrics = dict(key=key, V0=V0, initial=initial_name, forcing=input_name,
                           max_storage=float(storage.max()),
                           max_phase_error_rad=float(np.max(np.abs(zz[0]))),
                           max_PLL_frequency_hz=float(np.max(np.abs(freq))),
                           max_filtered_PLL_rocof_hz_s=float(np.max(np.abs(rocof))))
            assert metrics["max_storage"] <= 1+2e-6
            assert metrics["max_PLL_frequency_hz"] <= p["frequency_max_hz"]+1e-6
            assert metrics["max_filtered_PLL_rocof_hz_s"] <= p["rocof_max_hz_s"]+1e-6
            trajectories.append(metrics)
            histories[key] = np.column_stack([t, storage, freq, rocof, zz[0]])
np.savez_compressed(OUT/"nonlinear_histories.npz", **histories)
result = dict(status="NONLINEAR_IMPLEMENTATION_CHECK_PASSED_NOT_THE_SOURCE_OF_CERTIFICATION",
              random_points=5000, max_dissipation_residual=float(max_residual),
              max_sampled_boundary_derivative=float(max_boundary_derivative),
              trajectories=trajectories,
              synthesis_sha256=hashlib.sha256((OUT/"synthesis.json").read_bytes()).hexdigest(),
              verification_sha256=hashlib.sha256((OUT/"verification.json").read_bytes()).hexdigest(),
              validation_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/"nonlinear_validation.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
