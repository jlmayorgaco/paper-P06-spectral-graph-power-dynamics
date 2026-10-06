"""Independent exact-rational verification of nonlinear PLL sector certificates.

Rebuilds matrices from scalar model parameters, not saved rounded LMI matrices.
The saved IEEE-754 parameters and certificate X entries are interpreted as exact
binary rational values. The sine sector is proved analytically in the report.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import numpy as np
from pll_sector_contract import Specification, solve_contract

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/pll_sector_contract"


def rational(a):
    return F.from_float(float(a))


def zeros(n, m):
    return [[F(0) for j in range(m)] for i in range(n)]


def transpose(a):
    return [list(r) for r in zip(*a)]


def multiply(a, b):
    return [[sum((aa*bb for aa, bb in zip(row, col)), F(0))
             for col in zip(*b)] for row in a]


def exact_ldl_pivots(a):
    n = len(a)
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    ll = zeros(n, n)
    dd = []
    for i in range(n):
        ll[i][i] = F(1)
        pivot = a[i][i]-sum((ll[i][k]**2*dd[k] for k in range(i)), F(0))
        assert pivot > 0, f"nonpositive exact LDL pivot {i}: {float(pivot)}"
        dd.append(pivot)
        for j in range(i+1, n):
            ll[j][i] = (a[j][i]-sum((ll[j][k]*ll[i][k]*dd[k]
                                     for k in range(i)), F(0)))/pivot
    return dd


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def dual_psd(a):
    "Repair and verify each dual matrix, without assuming floating-point PSD."
    aa = [[(rational(a[i][j])+rational(a[j][i]))/2
           for j in range(len(a))] for i in range(len(a))]
    shift = F(1, 10**14)
    for attempt in range(12):
        b = [[aa[i][j]+(shift if i == j else 0) for j in range(len(a))]
             for i in range(len(a))]
        try:
            exact_ldl_pivots(b)
            return b, shift
        except AssertionError:
            shift *= 10
    raise AssertionError("could not repair dual PSD matrix")


def atan_interval(x, terms):
    total = sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
    next_term = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(total, total+next_term), max(total, total+next_term)


def verify(candidate, spec):
    p = {k:rational(v) for k,v in spec.items()}
    # Keep the normalized coordinates used in synthesis, whose frequency scale
    # is conservatively below the exact 2*pi*fmax by the explicit pi lower bound.
    sf = rational(2*np.pi*spec["frequency_max_hz"])
    scale = [p["angle_max"], sf, sf, sf]
    pi_lower = F("3.14159265358979323846264338327950288")
    atan5, atan239 = atan_interval(F(1, 5), 40), atan_interval(F(1, 239), 12)
    assert pi_lower <= 16*atan5[0]-4*atan239[1]  # Machin identity, exact bounds.
    assert sf <= 2*pi_lower*p["frequency_max_hz"]
    kp, ki = rational(candidate["Kp"]), rational(candidate["Ki"])
    eps = rational(candidate["epsilon_sdp"])
    X = [[rational(v) for v in row] for row in candidate["X_scaled"]]
    xpivots = exact_ldl_pivots(X)
    aends = [p["voltage_min"]*(1-p["angle_max"]**2/6), p["voltage_max"]]
    alpha = p["decay_rate"]*p["time_scale"]
    vertex_pivots = []
    rebuilt = []
    for a in aends:
        A = [[0, 1, 0, 0], [-kp*a/p["pll_tau"], -1/p["pll_tau"], 1/p["pll_tau"], 0],
             [-ki*a, 0, 0, 0], [0, 1/p["measurement_tau"], 0, -1/p["measurement_tau"]]]
        b = [[0, -p["frame_rate_per_voltage"]], [kp/p["pll_tau"], 0], [ki, 0], [0, 0]]
        A = [[p["time_scale"]*A[i][j]*scale[j]/scale[i] for j in range(4)] for i in range(4)]
        b = [[p["time_scale"]*b[i][j]/scale[i] for j in range(2)] for i in range(4)]
        AX, XAT = multiply(A, X), multiply(X, transpose(A))
        negM = zeros(6, 6)
        for i in range(4):
            for j in range(4):
                negM[i][j] = -(AX[i][j]+XAT[i][j]+2*alpha*X[i][j])
            for j in range(2):
                negM[i][4+j] = negM[4+j][i] = -eps*b[i][j]
        for j in range(2):
            negM[4+j][4+j] = 2*alpha
        vertex_pivots.append(exact_ldl_pivots(negM))
        exact_ldl_pivots([[negM[i][j]-(p["matrix_guard"] if i == j else 0)
                           for j in range(6)] for i in range(6)])
        rebuilt.append((A, b))
    r = p["frequency_max_hz"]/(p["measurement_tau"]*p["rocof_max_hz_s"])
    outputs = [[1, 0, 0, 0], [0, 1, 0, 0], [0, r, 0, -r]]
    margins = [1-multiply(multiply([c], X), transpose([c]))[0][0] for c in outputs]
    assert all(v > 0 for v in margins)
    guard = p["matrix_guard"]
    exact_ldl_pivots([[X[i][j]-(guard if i == j else 0)
                       for j in range(4)] for i in range(4)])
    assert all(v >= guard for v in margins)
    assert trace(X) <= p["certificate_trace_bound"]
    assert eps <= p["input_radius_upper_bound"]

    # A rigorous upper bound on the ideal guarded SDP optimum, using repaired
    # PSD duals and exact residual bounds on the compact primal search domain.
    ZX, shiftX = dual_psd(candidate["dual_X_lower"])
    duals = [dual_psd(a) for a in candidate["dual_vertices"]]
    lam = [max(F(0), rational(a)) for a in candidate["dual_outputs"]]
    lamT = max(F(0), rational(candidate["dual_trace"]))
    lamE = max(F(0), rational(candidate["dual_epsilon_upper"]))
    residualX = [[ZX[i][j]-sum((lam[k]*outputs[k][i]*outputs[k][j]
                              for k in range(3)), F(0))-(lamT if i == j else 0)
                  for j in range(4)] for i in range(4)]
    residualE = 1-lamE
    constant = (-guard*trace(ZX)+(1-guard)*sum(lam)
                +lamT*p["certificate_trace_bound"]+lamE*p["input_radius_upper_bound"])
    for (A, b), (Z, shift) in zip(rebuilt, duals):
        Z11 = [row[:4] for row in Z[:4]]
        za, atz = multiply(Z11, A), multiply(transpose(A), Z11)
        for i in range(4):
            for j in range(4):
                residualX[i][j] -= za[i][j]+atz[i][j]+2*alpha*Z11[i][j]
        residualE -= 2*sum((Z[i][4+j]*b[i][j] for i in range(4) for j in range(2)), F(0))
        constant += 2*alpha*(Z[4][4]+Z[5][5])-guard*trace(Z)
    rho = max(sum((abs(v) for v in row), F(0)) for row in residualX)
    upper = constant+p["certificate_trace_bound"]*rho+p["input_radius_upper_bound"]*max(F(0), residualE)
    assert upper >= eps
    return dict(status="EXACT_RATIONAL_LMI_AND_OUTPUT_CHECK_PASSED",
                epsilon_verified=float(eps),
                epsilon_verified_exact=str(eps),
                reference_frequency_envelope_hz=float(eps*p["frame_rate_per_voltage"]/(2*pi_lower)),
                min_X_LDL_pivot=float(min(xpivots)),
                min_negative_vertex_LDL_pivot=float(min(map(min, vertex_pivots))),
                exact_positive_output_margins=[str(v) for v in margins],
                output_margins_float=[float(v) for v in margins],
                gain_pair=[candidate["Kp"], candidate["Ki"]],
                fixed_gain_guarded_SDP_upper_bound=float(upper),
                fixed_gain_guarded_SDP_upper_bound_exact=str(upper),
                fixed_gain_guarded_SDP_relative_gap=float((upper-eps)/upper),
                dual_psd_diagonal_repairs=[float(shiftX)]+[float(s) for _, s in duals],
                dual_stationarity_X_row_sum_bound=float(rho),
                dual_stationarity_epsilon_residual=float(residualE),
                input_envelope="||Delta u_dq||^2/epsilon^2 + nu_reference^2/(frame_rate_per_voltage*epsilon)^2 <= 1",
                meaning="All measurable input signals within this moving-frame ellipsoid; initial state in certified ellipsoid; PLL and declared instrument only")


def main():
    path = OUT/"synthesis.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    results = [verify(c, data["specification"]) for c in data["tests"]+[data["best"]]]
    spec = Specification(**data["specification"])
    c = data["tests"][0]
    grad_audit = []
    for j, key in enumerate(("Kp", "Ki")):
        pp = np.array([c["Kp"], c["Ki"]])
        h = 1e-4*pp[j]
        plus, minus = pp.copy(), pp.copy()
        plus[j] += h; minus[j] -= h
        rp, rm = solve_contract(spec, *plus), solve_contract(spec, *minus)
        fd = (rp["epsilon_sdp"]-rm["epsilon_sdp"])/(2*h)
        analytic = c["gradient_epsilon_gains"][j]
        relative = abs(fd-analytic)/max(abs(fd), abs(analytic), 1e-12)
        assert relative < 2e-3, (key, fd, analytic, relative)
        grad_audit.append(dict(parameter=key, dual_gradient=analytic,
                               central_difference=fd, relative_error=relative))
    out = dict(status="CONDITIONAL_NONLINEAR_PLL_CONTRACT_VERIFIED",
               certificates=results, gradient_audit=grad_audit,
               synthesis_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
               verification_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               scope="Not a certificate for the GFL current/DC subsystem, algebraic network, SG replacement fraction, or raw bus RoCoF")
    (OUT/"verification.json").write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
