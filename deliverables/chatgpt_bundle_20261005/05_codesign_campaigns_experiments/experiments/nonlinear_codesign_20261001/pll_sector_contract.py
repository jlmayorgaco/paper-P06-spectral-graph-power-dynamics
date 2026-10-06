"""Exact sine-sector, bounded-voltage contract for the installed third-order PLL.

This certifies the PLL and a declared first-order RoCoF instrument conditional
on a terminal-voltage input envelope. It does NOT certify the DC/current loop,
network voltages, bus-frequency measurements, or SG-to-GFL replacement capacity.
"""
from pathlib import Path
from dataclasses import dataclass, asdict
import argparse
import hashlib
import json
import warnings
import numpy as np
import cvxpy as cp
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/pll_sector_contract"


@dataclass(frozen=True)
class Specification:
    voltage_min: float = 0.95
    voltage_max: float = 1.10
    angle_max: float = 1.0
    frequency_max_hz: float = 0.5
    rocof_max_hz_s: float = 0.5
    measurement_tau: float = 0.1
    pll_tau: float = 1 / (300 * 2 * np.pi)
    decay_rate: float = 5.0
    time_scale: float = 0.1
    matrix_guard: float = 1e-7
    frame_rate_per_voltage: float = 2*np.pi*10
    certificate_trace_bound: float = 1.0
    input_radius_upper_bound: float = 1.0

    def validate(self):
        assert 0 < self.voltage_min <= self.voltage_max
        assert 0 < self.angle_max < np.sqrt(6)
        assert min(self.frequency_max_hz, self.rocof_max_hz_s,
                   self.measurement_tau, self.pll_tau, self.decay_rate,
                   self.time_scale, self.matrix_guard) > 0
        assert self.frame_rate_per_voltage >= 0
        assert self.certificate_trace_bound > 0 and self.input_radius_upper_bound > 0


def state_scale(spec):
    return np.diag([spec.angle_max] + [2*np.pi*spec.frequency_max_hz]*3)


def vertices(spec, kp, ki):
    scale = state_scale(spec)
    inverse = np.diag(1/np.diag(scale))
    # Rigorous global sine inequality, not an unbounded Taylor truncation:
    # 1-delta^2/6 <= sin(delta)/delta <= 1 for |delta|<=angle_max<sqrt(6).
    amin = spec.voltage_min*(1-spec.angle_max**2/6)
    data = []
    for a in (amin, spec.voltage_max):
        A = np.array([[0, 1, 0, 0],
                      [-kp*a/spec.pll_tau, -1/spec.pll_tau, 1/spec.pll_tau, 0],
                      [-ki*a, 0, 0, 0],
                      [0, 1/spec.measurement_tau, 0, -1/spec.measurement_tau]])
        # theta_error_dot = omega_PLL - omega_reference. Keeping this input
        # preserves physical common frequency under a moving angular reference.
        b = np.array([[0, -spec.frame_rate_per_voltage],
                      [kp/spec.pll_tau, 0], [ki, 0], [0, 0]])
        Akp = np.zeros((4, 4)); Akp[1, 0] = -a/spec.pll_tau
        Aki = np.zeros((4, 4)); Aki[2, 0] = -a
        bkp = np.zeros((4, 2)); bkp[1, 0] = 1/spec.pll_tau
        bki = np.zeros((4, 2)); bki[2, 0] = 1
        data.append(dict(a=a, A=spec.time_scale*inverse@A@scale,
                         b=spec.time_scale*inverse@b,
                         Akp=spec.time_scale*inverse@Akp@scale,
                         Aki=spec.time_scale*inverse@Aki@scale,
                         bkp=spec.time_scale*inverse@bkp,
                         bki=spec.time_scale*inverse@bki))
    return data


def output_rows(spec):
    rscale = spec.frequency_max_hz/(spec.measurement_tau*spec.rocof_max_hz_s)
    return np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, rscale, 0, -rscale]])


def solve_contract(spec, kp, ki):
    spec.validate()
    X = cp.Variable((4, 4), symmetric=True)
    epsilon = cp.Variable(nonneg=True)
    alpha = spec.decay_rate*spec.time_scale
    vs = vertices(spec, kp, ki)
    constraints = [X >> spec.matrix_guard*np.eye(4)]
    constraints += [cp.quad_form(c, X) <= 1-spec.matrix_guard
                    for c in output_rows(spec)]
    constraints += [cp.trace(X) <= spec.certificate_trace_bound,
                    epsilon <= spec.input_radius_upper_bound]
    matrix_constraints = []
    for v in vs:
        upper = v["A"]@X + X@v["A"].T + 2*alpha*X
        M = cp.bmat([[upper, epsilon*v["b"]],
                     [epsilon*v["b"].T, -2*alpha*np.eye(2)]])
        cc = M << -spec.matrix_guard*np.eye(6)
        constraints.append(cc); matrix_constraints.append(cc)
    problem = cp.Problem(cp.Maximize(epsilon), constraints)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            problem.solve(solver="CLARABEL", tol_gap_abs=1e-10,
                          tol_gap_rel=1e-10, tol_feas=1e-10,
                          max_iter=300, max_threads=1)
    except cp.error.SolverError as exc:
        return dict(status="SOLVER_ERROR", error=str(exc), Kp=kp, Ki=ki)
    if X.value is None or epsilon.value is None:
        return dict(status=problem.status, Kp=kp, Ki=ki)
    xx = (X.value+X.value.T)/2
    eps = float(epsilon.value)
    derivative = []
    matrices = []
    for v, constraint in zip(vs, matrix_constraints):
        matrices.append(np.block([
            [v["A"]@xx+xx@v["A"].T+2*alpha*xx, eps*v["b"]],
            [eps*v["b"].T, -2*alpha*np.eye(2)]]))
    for key in ("kp", "ki"):
        derivative.append(-sum(float(np.sum(cc.dual_value*np.block([
            [v["A"+key]@xx+xx@v["A"+key].T, eps*v["b"+key]],
            [eps*v["b"+key].T, np.zeros((2, 2))]])))
            for v, cc in zip(vs, matrix_constraints)))
    return dict(status=problem.status, Kp=float(kp), Ki=float(ki),
                epsilon_sdp=eps, X_scaled=xx.tolist(),
                reference_frequency_envelope_hz=eps*spec.frame_rate_per_voltage/(2*np.pi),
                gradient_epsilon_gains=derivative,
                max_vertex_eigenvalue=max(float(np.linalg.eigvalsh(m)[-1]) for m in matrices),
                min_X_eigenvalue=float(np.linalg.eigvalsh(xx)[0]),
                normalized_output_bounds=np.sqrt(np.einsum(
                    "ij,jk,ik->i", output_rows(spec), xx, output_rows(spec))).tolist(),
                dual_X_lower=constraints[0].dual_value.tolist(),
                dual_outputs=[float(c.dual_value) for c in constraints[1:4]],
                dual_trace=float(constraints[4].dual_value),
                dual_epsilon_upper=float(constraints[5].dual_value),
                dual_vertices=[c.dual_value.tolist() for c in matrix_constraints],
                solver="CLARABEL", solve_time=problem.solver_stats.solve_time)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--optimize", action="store_true")
    args = parser.parse_args()
    spec = Specification()
    OUT.mkdir(parents=True, exist_ok=True)
    nominal = np.array([2*np.pi*5, (2*np.pi*5)**2/4])
    tests = [solve_contract(spec, *nominal), solve_contract(spec, 40.0, 355.438)]
    history = []
    optimization = None
    if args.optimize:
        cache = {}

        def oracle(loggains):
            key = tuple(loggains)
            if key not in cache:
                r = solve_contract(spec, *np.exp(loggains))
                cache[key] = r
                history.append(r)
                print(json.dumps({k: r.get(k) for k in ("status", "Kp", "Ki", "epsilon_sdp")}), flush=True)
            r = cache[key]
            if "epsilon_sdp" not in r:
                return 1.0, np.zeros(2)
            return -r["epsilon_sdp"], -np.array(r["gradient_epsilon_gains"])*np.exp(loggains)

        opt = minimize(oracle, np.log(nominal), jac=True, method="L-BFGS-B",
                       bounds=list(zip(np.log(0.25*nominal), np.log(4*nominal))),
                       options=dict(maxiter=45, ftol=1e-10, gtol=1e-7, maxls=20))
        optimization = dict(success=bool(opt.success), message=str(opt.message),
                            iterations=opt.nit, evaluations=opt.nfev,
                            gains=np.exp(opt.x).tolist(),
                            status="LOCAL_GAIN_SEARCH_NOT_GLOBAL_OPTIMUM")
    candidates = [r for r in tests+history if "epsilon_sdp" in r and
                  r["min_X_eigenvalue"] > 0 and r["max_vertex_eigenvalue"] < 1e-7]
    if not candidates:
        result = dict(status="NO_FEASIBLE_CONTRACT_FOUND", specification=asdict(spec), tests=tests)
    else:
        best = max(candidates, key=lambda r:r["epsilon_sdp"])
        result = dict(status="SDP_CANDIDATE_REQUIRES_INDEPENDENT_CERTIFICATE_CHECK",
                      specification=asdict(spec), tests=tests, best=best,
                      optimization=optimization, history=history,
                      cvxpy_version=cp.__version__, numpy_version=np.__version__,
                      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      scope="PLL plus filtered RoCoF instrument; joint moving-frame voltage/reference-rate ellipsoid, no network or full-converter certificate")
    (OUT/"synthesis.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: result.get(k) for k in ("status", "tests", "best", "optimization")}, indent=2))


if __name__ == "__main__":
    main()
