"""Full installed GFL + measurement: exact nonlinear polytopic containment.

The state contains PLL, current-controller integrators, filter currents, DC
controller, DC voltage and a causal RoCoF instrument. Local inputs remain ports
to be closed by the network; no replacement-capacity certificate is asserted.
"""
from pathlib import Path
from dataclasses import dataclass, asdict
from itertools import product
import argparse
import hashlib
import json
import time
import tomllib
import warnings
import numpy as np
import cvxpy as cp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"reports/nonlinear_codesign_20261001/full_gfl_sector_contract"
STATE_NAMES = ["delta", "omega", "xi", "chi", "gamma_q", "gamma_d",
               "i_d", "i_q", "v_dc_i", "v_dc"]


@dataclass(frozen=True)
class Specification:
    angle_box: float = 0.05
    id_box: float = 0.1
    iq_box: float = 0.1
    vdc_box: float = 0.05
    gamma_scale: float = 0.001
    integral_current_scale: float = 0.1
    frequency_max_hz: float = 0.5
    rocof_max_hz_s: float = 0.5
    measurement_tau: float = 0.1
    decay_rate: float = 5.0
    time_scale: float = 0.1
    frame_rate_per_voltage: float = 2*np.pi*10
    dc_power_per_voltage: float = 10.0
    matrix_guard: float = 1e-8
    trace_bound: float = 2.0
    epsilon_bound: float = 1.0
    objective_radius_unit: float = 0.001


def state_scale(s):
    return np.array([s.angle_box]+[2*np.pi*s.frequency_max_hz]*3+
                    [s.gamma_scale]*2+[s.id_box, s.iq_box,
                    s.integral_current_scale, s.vdc_box])


def parameter_bounds(s, d):
    assert 0 < s.angle_box < np.sqrt(6)
    Vdc = d["parameters"]["Vdc"]
    assert 0 < s.vdc_box < Vdc
    return [(1-s.angle_box**2/6, 1.0), (-s.angle_box/2, s.angle_box/2),
            (-s.id_box, s.id_box), (-s.iq_box, s.iq_box),
            (1/(Vdc+s.vdc_box), 1/(Vdc-s.vdc_box))]


def physical_matrices(s, d, kp, ki, parameters):
    "Return exact factorization A(z) z + epsilon B(z) eta; no Jacobian truncation."
    ss, cc, did, diq, beta = parameters
    p = d["parameters"]
    A = np.zeros((10, 10)); B = np.zeros((10, 4))
    a, omega0 = p["omega_base"]/p["Xf"], p["omega_base"]*p["omega_frame"]
    cd, cq = np.zeros(10), np.zeros(10)
    cd[5], cd[6], cd[8], cd[9] = p["cc_ki"], -p["cc_kp"], p["cc_kp"], p["cc_kp"]*p["dc_kp"]
    cq[4], cq[7] = p["cc_ki"], -p["cc_kp"]
    A[0, 1] = 1; B[0, 2] = -s.frame_rate_per_voltage
    A[1, 0], A[1, 1], A[1, 2] = -kp*d["V0"]*ss/p["pll_tau"], -1/p["pll_tau"], 1/p["pll_tau"]
    B[1, 1] = kp/p["pll_tau"]
    A[2, 0] = -ki*d["V0"]*ss; B[2, 1] = ki
    A[3, 1], A[3, 3] = 1/s.measurement_tau, -1/s.measurement_tau
    A[4, 7] = -1
    A[5, 6], A[5, 8], A[5, 9] = -1, 1, p["dc_kp"]
    A[6] = a*cd; A[6, 6] -= a*p["Rf"]; A[6, 7] += omega0
    A[6, 0] -= a*d["V0"]*cc
    A[6, 1] += d["iq0"]+diq; B[6, 0] = -a
    A[7] = a*cq; A[7, 7] -= a*p["Rf"]; A[7, 6] -= omega0
    A[7, 0] += a*d["V0"]*ss
    A[7, 1] -= d["id0"]+did; B[7, 1] = -a
    A[8, 9] = p["dc_ki"]
    hh = d["id0"]*cd+d["iq0"]*cq
    hh[6] += d["vid0"]; hh[7] += d["viq0"]
    A[9] = -beta/p["Cdc"]*(hh+did*cd+diq*cq)
    B[9, 3] = beta/p["Cdc"]*s.dc_power_per_voltage
    return A, B


def exact_rhs(s, d, kp, ki, z, epsilon, eta):
    "Direct nonlinear RHS in PLL-aligned current coordinates for independent checks."
    delta, omega, xi, chi, gq, gd, did, diq, vdi, dvdc = z
    p = d["parameters"]
    ud = d["V0"]*np.cos(delta)+epsilon*eta[0]
    uq = -d["V0"]*np.sin(delta)+epsilon*eta[1]
    id_, iq = d["id0"]+did, d["iq0"]+diq
    ed, eq = p["dc_kp"]*dvdc+vdi-did, -diq
    vid = d["vid0"]+p["cc_kp"]*ed+p["cc_ki"]*gd
    viq = d["viq0"]+p["cc_kp"]*eq+p["cc_ki"]*gq
    freq = p["omega_base"]*p["omega_frame"]+omega
    a = p["omega_base"]/p["Xf"]
    # Incremental power expression avoids using an approximately trimmed constant.
    dpac = (vid-d["vid0"])*d["id0"]+(viq-d["viq0"])*d["iq0"]
    dpac += vid*did+viq*diq
    return np.array([omega-epsilon*s.frame_rate_per_voltage*eta[2],
        (xi+kp*uq-omega)/p["pll_tau"], ki*uq, (omega-chi)/s.measurement_tau,
        eq, ed,
        a*((vid-d["vid0"])-(ud-d["V0"])-p["Rf"]*did)
             +p["omega_base"]*p["omega_frame"]*diq+omega*iq,
        a*((viq-d["viq0"])-uq-p["Rf"]*diq)
             -p["omega_base"]*p["omega_frame"]*did-omega*id_,
        p["dc_ki"]*dvdc,
        (epsilon*s.dc_power_per_voltage*eta[3]-dpac)/(p["Cdc"]*(p["Vdc"]+dvdc))])


def state_parameters(s, d, z):
    delta = z[0]
    ss = np.sinc(delta/np.pi)
    # Stable evaluation avoids cancellation in (cos(delta)-1)/delta.
    cc = -2*np.sin(delta/2)**2/delta if delta != 0 else 0.0
    return [ss, cc, z[6], z[7], 1/(d["parameters"]["Vdc"]+z[9])]


def output_rows(s):
    C = np.eye(10)[[0, 1, 6, 7, 9], :]
    cr = np.zeros(10)
    cr[1] = s.frequency_max_hz/(s.measurement_tau*s.rocof_max_hz_s)
    cr[3] = -cr[1]
    return np.vstack([C, cr])


def vertices(s, d, kp, ki):
    scale = state_scale(s)
    result = []
    for pars in product(*parameter_bounds(s, d)):
        A, B = physical_matrices(s, d, kp, ki, pars)
        result.append((s.time_scale*A*scale[None, :]/scale[:, None],
                       s.time_scale*B/scale[:, None], list(pars)))
    return result


def solve_contract(s, d, kp, ki, formulation="radius"):
    tstart = time.monotonic()
    vs = vertices(s, d, kp, ki)
    alpha = s.decay_rate*s.time_scale
    worst_abscissa = max(float(np.linalg.eigvals(A).real.max()/s.time_scale) for A, B, pp in vs)
    if worst_abscissa >= -s.decay_rate:
        return dict(status="SECTOR_VERTEX_DECAY_INFEASIBLE", bus=d["bus"],
                    Kp=kp, Ki=ki, worst_vertex_abscissa=worst_abscissa)
    X = cp.Variable((10, 10), symmetric=True)
    radius = cp.Variable(nonneg=True)
    squared = formulation == "squared_radius"
    eps = s.objective_radius_unit*radius if not squared else None
    radius_squared = s.objective_radius_unit**2*radius if squared else None
    C = output_rows(s)
    constraints = [X >> s.matrix_guard*np.eye(10)]
    constraints += [cp.quad_form(c, X) <= 1-s.matrix_guard for c in C]
    constraints += [cp.trace(X) <= s.trace_bound,
                    radius_squared <= s.epsilon_bound**2 if squared else eps <= s.epsilon_bound]
    matrix_constraints = []
    for A, B, pp in vs:
        if squared:
            # Exact Schur complement of the guarded 14x14 inequality. This
            # optimizes epsilon^2 in a 10x10 affine LMI and improves scaling.
            N = (A@X+X@A.T+2*alpha*X+s.matrix_guard*np.eye(10)
                 +radius_squared/(2*alpha-s.matrix_guard)*(B@B.T))
            con = N << 0
        else:
            M = cp.bmat([[A@X+X@A.T+2*alpha*X, eps*B], [eps*B.T, -2*alpha*np.eye(4)]])
            con = M << -s.matrix_guard*np.eye(14)
        constraints.append(con); matrix_constraints.append(con)
    problem = cp.Problem(cp.Maximize(radius), constraints)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            problem.solve(solver="CLARABEL", tol_gap_abs=2e-10, tol_gap_rel=2e-10,
                          tol_feas=2e-10, max_iter=250, max_threads=1)
    except cp.error.SolverError as exc:
        return dict(status="SOLVER_ERROR", error=str(exc), bus=d["bus"], Kp=kp, Ki=ki)
    if X.value is None or radius.value is None:
        return dict(status=problem.status, bus=d["bus"], Kp=kp, Ki=ki,
                    worst_vertex_abscissa=worst_abscissa)
    xx = (X.value+X.value.T)/2
    epsilon = s.objective_radius_unit*np.sqrt(float(radius.value)) if squared else float(eps.value)
    dual_scale = s.objective_radius_unit**2/(2*epsilon) if squared else s.objective_radius_unit
    dual_vertices=[]
    for (A,B,pp),con in zip(vs,matrix_constraints):
        if squared:
            lift=np.vstack([np.eye(10),epsilon/(2*alpha-s.matrix_guard)*B.T])
            dual_vertices.append((dual_scale*lift@con.dual_value@lift.T).tolist())
        else:
            dual_vertices.append((dual_scale*con.dual_value).tolist())
    matrices = [np.block([[A@xx+xx@A.T+2*alpha*xx, epsilon*B],
                          [epsilon*B.T, -2*alpha*np.eye(4)]]) for A, B, pp in vs]
    return dict(status=problem.status, bus=d["bus"], Kp=float(kp), Ki=float(ki),
        epsilon_sdp=epsilon, X_scaled=xx.tolist(),
        max_vertex_eigenvalue=max(float(np.linalg.eigvalsh(M)[-1]) for M in matrices),
        min_X_eigenvalue=float(np.linalg.eigvalsh(xx)[0]),
        normalized_output_bounds=np.sqrt(np.einsum("ij,jk,ik->i",C,xx,C)).tolist(),
        worst_vertex_abscissa=worst_abscissa,
        dual_X_lower=(dual_scale*constraints[0].dual_value).tolist(),
        dual_outputs=[dual_scale*float(c.dual_value) for c in constraints[1:1+len(C)]],
        dual_trace=dual_scale*float(constraints[1+len(C)].dual_value),
        dual_epsilon_upper=dual_scale*float(constraints[2+len(C)].dual_value)*(2*s.epsilon_bound if squared else 1),
        dual_vertices=dual_vertices, formulation=formulation,
        solve_time=problem.solver_stats.solve_time, elapsed=time.monotonic()-tstart)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buses", default="all")
    parser.add_argument("--gain-policy", choices=["nominal", "network_candidate", "pll_contract"], default="network_candidate")
    parser.add_argument("--label", default="initial")
    parser.add_argument("--formulation",choices=["squared_radius","radius"],default="radius")
    args = parser.parse_args()
    path = OUT/"models.toml"
    models = tomllib.loads(path.read_text(encoding="utf-8"))
    source = tomllib.loads((OUT.parent/"candidate_final_physical.toml").read_text(encoding="utf-8"))
    selected = None if args.buses == "all" else [int(b) for b in args.buses.split(",")]
    s = Specification()
    result = dict(status="FULL_GFL_SDP_CANDIDATES_REQUIRE_INDEPENDENT_VERIFICATION", formulation=args.formulation,
                  specification=asdict(s), state_names=STATE_NAMES, gain_policy=args.gain_policy,
                  models_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  candidates=[], scope="Full local GFL plus filtered PLL RoCoF; network ports remain open")
    for d in models["devices"]:
        if selected is not None and d["bus"] not in selected:
            continue
        if args.gain_policy == "nominal":
            kp, ki = d["Kp_nominal"], d["Ki_nominal"]
        elif args.gain_policy == "network_candidate":
            idx = d["bus"]-30; kp, ki = source["Kp"][idx], source["Ki"][idx]
        else:
            pp = json.loads((OUT.parent/"pll_sector_contract/synthesis.json").read_text(encoding="utf-8"))["best"]
            kp, ki = pp["Kp"], pp["Ki"]
        r = solve_contract(s, d, kp, ki,args.formulation)
        result["candidates"].append(r)
        print(json.dumps({k:r.get(k) for k in ("bus", "status", "epsilon_sdp", "worst_vertex_abscissa", "elapsed")}), flush=True)
        (OUT/f"synthesis_{args.label}.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
