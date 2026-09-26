"""Journal gate 2 - nonlinear phasor-domain time-domain validation.

This is NOT an electromagnetic-transient (EMT) study. It is a nonlinear
phasor-domain DAE simulation: the network is its algebraic current balance,
the devices are their differential equations, integrated without linearization
(BDF, algebraic equations solved by Newton at every right-hand-side call, the E25
and E32 integrator). It tests whether the modes and the stability verdicts of the
eigenanalysis are what the nonlinear model actually does.

Everything below was declared before the full run. Cases, horizons and outcome
rules have not changed since the first test runs. The disturbance and the
integrator's Jacobian were changed after test runs, for the reasons recorded
below.

Cases (members, operating point, prediction from the linear model):

    IEEE-39 (F7A slice, t = 1.5, h = 1; leaky Q/V regulator)
      SAFE        P_inf flagship 30+33+35+37      g = 1,       k = 0.5
      K4          P4 flagship                      g = 0.03625, k = 1.425
      K4-proper   P4 triple 30+33+35               (stable when kappa = 4)
      K2          P2 pair 30+33                    g = 0.11,    k = 1.85
      K2-single   P2 singles 30 and 33             (stable when kappa = 2)
      TONGUE      flagship on the line k = 1.30:   g = 0.02 (before the gap),
                  0.042 (gap), 0.102 (tongue interior), 0.173 (after the tongue)
      COND-       P4 flagship + condenser D = 2 (RHP-validated, G1) at 2.0 % rating
      COND+       the same at 3.0 % (RHP threshold 2.48 %)
      SWING       P_inf, replacement 37 + undamped classical condenser at 25 %
                  (the out-of-band mode found by G1)
    Kundur (K12A slice, t = 1)
      OSC-        replacement 2, k = 1.25, g = 0.08   (RHP boundary at g = 0.0949)
      OSC+        replacement 2, k = 1.25, g = 0.11
      APER        replacements 2+3, k = 0.75, g = 0    (real eigenvalue)

Disturbances (identical for every case of a grid; both leave the power balance,
and hence the equilibrium, unchanged):
    D1  rotor-speed kick +1e-4 pu on a machine present in every case
        (IEEE-39: bus 31; Kundur: machine 1)
    D2  +2 % active-load pulse lasting 0.2 s (IEEE-39: bus 20; Kundur: bus 7),
        integrated as two segments so that the solver never steps across the switch
A sustained load step (the E25/E32 disturbance) was tried first and rejected in
test runs, before any case was run in full. Without governors a sustained step
has no equilibrium: the common frequency drifts without bound (-0.008 pu after
10 s on Kundur), leaves the small-signal regime and ends in a numerical
collapse. The linear analysis excludes that drift by removing the reference mode
(FAILED_EXPERIMENTS F18).

Horizon T = clip(3 / |alpha_c|, 30, 120) s, with alpha_c the rightmost
non-reference eigenvalue. The run stops early if the relative machine-speed
spread exceeds 0.05 pu (DIVERGED).

Observable: the difference of the two machine speeds whose components in the
right eigenvector of the predicted critical mode differ most (a model-informed
choice of WHICH physical signal to look at, not of what it does). For a real
critical eigenvalue, the deviation of the state with the largest participation
instead, with the growth rate fitted to its logarithm until it reaches 0.05; a
real mode need not move machines against each other.

Integrator: BDF, with a central-difference Jacobian of the reduced right-hand
side (fixed relative step 1e-6). Two earlier choices were replaced:
- The equilibrium A_red as a fixed Jacobian (E25/E32). In Kundur test runs,
  BDF's Newton iteration stalled on fast converter PLL transients.
- scipy's adaptive num_jac, used in the first full run. It grew its step on
  columns that are identically zero (the frozen EMFs of the condenser cases),
  hit states with no network solution, and returned NaN. The resulting LU error
  was misreported as ALGEBRAIC_FAILURE in the condenser and swing cases. That
  run is kept as results/G2/first_run/.
Every case was re-run with the final integrator. Cases, disturbances, horizons
and outcome rules are unchanged.

Model-free outcome from the trajectory:
    DIVERGED   spread > 0.05 pu, or the algebraic equations have no solution
    UNSTABLE   envelope growth rate of the second half > 0
    STABLE     envelope growth rate of the second half < 0
Frequency: FFT peak of the second half (0.05-20 Hz), and the dominant
matrix-pencil exponent (largest contribution at the end of the small-signal
window, spread <= 1e-3 pu); growth: the same pencil exponent and the envelope.
"""

from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.integrate import BDF, OdeSolution
from scipy.linalg import lu_factor, lu_solve
from scipy.signal import find_peaks

from _bootstrap import RESULTS
from _f7_common import Theta
from _gates import GateExperiment
from _overnight import pin_blas_threads
from E25_time_domain_validation import matrix_pencil
from F8_service_attribution import r_configs, solve_config
from F12_kundur import solve as k_solve
from G1_f8_rhp_followup import _em_config
from ibr_cycles.dynamics.linearize import central_difference_jacobians

OUT = RESULTS / "G2"
STEP = 0.02
DIVERGE = 0.05
LINEAR = 1e-3
T_MIN, T_MAX = 30.0, 120.0
FIT_START = 5.0
F7A = {"t": 1.5, "h": 1.0}
FLAG = (30, 33, 35, 37)
KICK = 1e-4
PULSE_S = 0.2
DISTURBANCES = {
    "IEEE-39": {
        "D1": {"kind": "kick", "machine": 31, "load_bus": 20},
        "D2": {"kind": "pulse", "machine": None, "load_bus": 20},
    },
    "Kundur": {
        "D1": {"kind": "kick", "machine": 1, "load_bus": 7},
        "D2": {"kind": "pulse", "machine": None, "load_bus": 7},
    },
}

CASES = [
    ("IEEE-39", "SAFE P_inf flagship", FLAG, {"g": 1.0, "k": 0.5}, None),
    ("IEEE-39", "K4 P4 flagship", FLAG, {"g": 0.03625, "k": 1.425}, None),
    (
        "IEEE-39",
        "K4-proper P4 30+33+35",
        (30, 33, 35),
        {"g": 0.03625, "k": 1.425},
        None,
    ),
    ("IEEE-39", "K2 P2 30+33", (30, 33), {"g": 0.11, "k": 1.85}, None),
    ("IEEE-39", "K2-single P2 30", (30,), {"g": 0.11, "k": 1.85}, None),
    ("IEEE-39", "K2-single P2 33", (33,), {"g": 0.11, "k": 1.85}, None),
    ("IEEE-39", "TONGUE g=0.020 before gap", FLAG, {"g": 0.020, "k": 1.30}, None),
    ("IEEE-39", "TONGUE g=0.042 gap", FLAG, {"g": 0.042, "k": 1.30}, None),
    ("IEEE-39", "TONGUE g=0.102 interior", FLAG, {"g": 0.102, "k": 1.30}, None),
    ("IEEE-39", "TONGUE g=0.173 after", FLAG, {"g": 0.173, "k": 1.30}, None),
    (
        "IEEE-39",
        "COND- P4 flagship + SC(D=2) 2.0%",
        FLAG,
        {"g": 0.03625, "k": 1.425},
        ("R_0010000", 0.020),
    ),
    (
        "IEEE-39",
        "COND+ P4 flagship + SC(D=2) 3.0%",
        FLAG,
        {"g": 0.03625, "k": 1.425},
        ("R_0010000", 0.030),
    ),
    (
        "IEEE-39",
        "SWING P_inf 37 + undamped classical SC 25%",
        (37,),
        {"g": 1.0, "k": 0.5},
        ("R_0000000", 0.25),
    ),
    ("Kundur", "OSC- K12A 2 g=0.08", (2,), {"g": 0.08, "k": 1.25}, None),
    ("Kundur", "OSC+ K12A 2 g=0.11", (2,), {"g": 0.11, "k": 1.25}, None),
    ("Kundur", "APER K12A 2+3 g=0", (2, 3), {"g": 0.0, "k": 0.75}, None),
]


def build(bench, members, point, condenser):
    if bench == "Kundur":
        return k_solve(members, {**point, "t": 1.0})
    theta = Theta(g=point["g"], k=point["k"], **F7A)
    configs = {c["name"]: c for c in r_configs()}
    if condenser is None:
        return solve_config(members, theta, configs["R_none"])
    name, rating = condenser
    return solve_config(members, theta, _em_config(configs[name], rating))


def critical(case):
    values, right = np.linalg.eig(case.system.A)
    keep = np.abs(values) > 1e-3
    idx = np.flatnonzero(keep & (values.imag >= 0))
    c = idx[np.argmax(values[idx].real)]
    rhp = int(np.count_nonzero(keep & (values.real > 0) & (values.imag >= 0)))
    return values[c], right[:, c], rhp


def observable(case, vector):
    labels = case.system.labels
    speeds = [i for i, n in enumerate(labels) if n.startswith("omega_sg")]
    comp = vector[speeds]
    best, pair = -1.0, (speeds[0], speeds[-1])
    for a in range(len(speeds)):
        for b in range(a + 1, len(speeds)):
            d = abs(comp[a] - comp[b])
            if d > best:
                best, pair = d, (speeds[a], speeds[b])
    return pair, speeds


class AlgebraicFailure(RuntimeError):
    """The network equations have no nearby solution (voltage collapse)."""


def simulate(case, bus, horizon, speeds):
    """BDF stepped by hand, so divergence or collapse keeps the trajectory so far."""

    dae = case.dae
    network = dae.network
    x0, z0 = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    state = {
        "z": z0.copy(),
        "factor": lu_factor(central_difference_jacobians(dae, x0, z0, {}).gz),
    }
    original = network.loads[bus["load_bus"]]
    if bus["kind"] == "kick":
        x0[case.system.labels.index(f"omega_sg{bus['machine']}")] += KICK

    def gz(x, z, h=1e-7):
        cols = []
        for k in range(z.size):
            dz = np.zeros_like(z)
            dz[k] = h
            cols.append((dae.g(x, z + dz, {}) - dae.g(x, z - dz, {})) / (2 * h))
        return np.column_stack(cols)

    def rhs(t, x):
        # Chord iterations from the last converged z; if they stall, full Newton
        # (gz rebuilt every iteration), whose last factorization becomes the new
        # chord matrix. A trial point where the network equations have no solution
        # returns NaN, which makes BDF shrink the step; only if no step size
        # succeeds does the run end (SOLVER_FAILED = collapse).
        z = state["z"].copy()
        for _ in range(15):
            residual = dae.g(x, z, {})
            if not np.all(np.isfinite(residual)):
                break
            if np.abs(residual).max() < 1e-10:
                state["z"] = z
                return dae.f(x, z, {})
            z = z - lu_solve(state["factor"], residual)
        z = state["z"].copy()
        for _ in range(12):
            residual = dae.g(x, z, {})
            if not np.all(np.isfinite(residual)):
                break
            if np.abs(residual).max() < 1e-10:
                state["z"] = z
                return dae.f(x, z, {})
            try:
                state["factor"] = lu_factor(gz(x, z))
            except (np.linalg.LinAlgError, ValueError):
                break
            z = z - lu_solve(state["factor"], residual)
        return np.full(x.shape, np.nan)

    def jac(t, x):
        # Central differences with a FIXED relative step. scipy's adaptive num_jac
        # grows its step on columns whose derivative is identically zero (frozen
        # condenser EMFs), reaches states where the network equations have no
        # solution, and returns NaN; that was misreported as a collapse.
        n = x.size
        out = np.empty((n, n))
        for j in range(n):
            h = 1e-6 * max(1.0, abs(x[j]))
            up, down = x.copy(), x.copy()
            up[j] += h
            down[j] -= h
            out[:, j] = (rhs(t, up) - rhs(t, down)) / (2.0 * h)
        return out if np.all(np.isfinite(out)) else case.system.A

    # a load pulse is two segments, so the solver never steps across the switch
    segments = [(0.0, horizon, False)]
    if bus["kind"] == "pulse":
        segments = [(0.0, PULSE_S, True), (PULSE_S, horizon, False)]
    times, pieces, status, x = [0.0], [], "COMPLETED", x0
    try:
        for t0, t1, stepped in segments:
            network.loads[bus["load_bus"]] = (
                complex(original.real * (1.0 + STEP), original.imag)
                if stepped
                else original
            )
            solver = BDF(rhs, t0, x, t1, jac=jac, rtol=1e-7, atol=1e-9)
            while solver.status == "running":
                try:
                    solver.step()
                except (AlgebraicFailure, np.linalg.LinAlgError, ValueError):
                    status = "ALGEBRAIC_FAILURE"
                    break
                if solver.status == "failed" or not np.all(np.isfinite(solver.y)):
                    status = "SOLVER_FAILED"
                    break
                times.append(solver.t)
                pieces.append(solver.dense_output())
                if np.abs(solver.y[speeds] - solver.y[speeds].mean()).max() > DIVERGE:
                    status = "DIVERGED"
                    break
            if status != "COMPLETED":
                break
            x = solver.y.copy()
    finally:
        network.loads[bus["load_bus"]] = original
    return OdeSolution(times, pieces), float(times[-1]), status


def dominant_exponent(times, signal):
    """Pencil exponents and least-squares amplitudes; the one dominating the end."""

    ex = matrix_pencil(times, signal, order=14)
    tau = times - times[0]
    # a component growing or decaying by more than e^30 over the window is a
    # numerical artefact of the pencil, not a mode of the record
    ex = ex[np.isfinite(ex) & (np.abs(ex.real) * tau[-1] < 30.0)]
    if ex.size == 0:
        return complex("nan"), ex
    y = signal - signal.mean()
    basis = np.exp(np.outer(tau, ex))
    amp, *_ = np.linalg.lstsq(basis, y.astype(complex), rcond=None)
    weight = np.abs(amp) * np.exp(ex.real * tau[-1])
    order = np.argsort(-weight)
    for i in order:
        if ex[i].imag >= 0:
            return complex(ex[i]), ex
    return complex(ex[order[0]]), ex


def envelope_growth(times, signal):
    y = np.abs(signal - signal.mean())
    peaks, _ = find_peaks(y)
    if peaks.size >= 4:
        return float(
            np.polyfit(times[peaks], np.log(np.maximum(y[peaks], 1e-300)), 1)[0]
        )
    return float(np.polyfit(times, np.log(np.maximum(y, 1e-300)), 1)[0])


def fft_peak(times, signal, lo=0.05, hi=20.0):
    y = (signal - signal.mean()) * np.hanning(signal.size)
    spec = np.abs(np.fft.rfft(y))
    freq = np.fft.rfftfreq(signal.size, d=float(times[1] - times[0]))
    keep = (freq >= lo) & (freq <= hi)
    return float(freq[keep][np.argmax(spec[keep])]) if keep.any() else float("nan")


def run(task):
    bench, name, members, point, condenser, dist = task
    case = build(bench, members, point, condenser)
    lam, vec, rhp = critical(case)
    (i, j), speeds = observable(case, vec)
    labels = case.system.labels
    horizon = float(np.clip(3.0 / max(abs(lam.real), 1e-9), T_MIN, T_MAX))
    started = time.time()
    dspec = DISTURBANCES[bench][dist]
    sol, t_end, status = simulate(case, dspec, horizon, speeds)
    row = {
        "benchmark": bench,
        "case": name,
        "members": "+".join(map(str, members)),
        "disturbance": dist,
        "disturbance_kind": dspec["kind"],
        **point,
        "condenser": f"{condenser[0]}@{condenser[1]:g}" if condenser else "",
        "lambda_re": float(lam.real),
        "lambda_f_hz": float(abs(lam.imag) / (2 * math.pi)),
        "rhp_count": rhp,
        "predicted": "UNSTABLE" if rhp else "STABLE",
        "observable": f"{labels[i]}-{labels[j]}",
        "horizon_s": horizon,
        "t_end_s": t_end,
        "cpu_s": time.time() - started,
        "run_status": status,
    }
    n = int(min(6000, max(1500, t_end * 50)))
    times = np.linspace(0.0, t_end, n)
    x = sol(times)
    s = x[i] - x[j]
    spread = np.abs(x[speeds] - x[speeds].mean(axis=0)).max(axis=0)
    diverged = status in ("DIVERGED", "ALGEBRAIC_FAILURE", "SOLVER_FAILED")
    half = times >= 0.5 * t_end
    row["env_growth_2nd_half"] = (
        envelope_growth(times[half], s[half]) if half.sum() > 50 else float("nan")
    )
    row["fft_f_hz"] = (
        fft_peak(times[half], s[half]) if half.sum() > 50 else float("nan")
    )
    row["max_spread_pu"] = float(spread.max())
    in_linear = np.maximum.accumulate(spread) <= LINEAR
    t_lin = float(times[in_linear][-1]) if in_linear.any() else 0.0
    lin = in_linear & (times >= min(FIT_START, t_lin / 3.0))
    if lin.sum() > 200:
        tt = np.linspace(times[lin][0], times[lin][-1], min(3000, int(lin.sum())))
        ss = np.interp(tt, times, s)
        dom, ex = dominant_exponent(tt, ss)
        closest = ex[np.argmin(np.abs(ex - lam))] if ex.size else complex("nan")
        row.update(
            pencil_window_s=float(tt[-1] - tt[0]),
            pencil_re=float(dom.real),
            pencil_f_hz=float(abs(dom.imag) / (2 * math.pi)),
            closest_re=float(closest.real),
            closest_f_hz=float(abs(closest.imag) / (2 * math.pi)),
        )
        # the rightmost inter-area (0.3-1.5 Hz) eigenvalue and its pencil match:
        # the mode that crosses at the condenser and policy thresholds, reported
        # even when another mode is rightmost
        allv = np.linalg.eigvals(case.system.A)
        fb = np.abs(allv.imag) / (2 * math.pi)
        band = allv[(fb >= 0.3) & (fb <= 1.5) & (allv.imag > 0)]
        if band.size and ex.size:
            lb = band[np.argmax(band.real)]
            mb = ex[np.argmin(np.abs(ex - lb))]
            row.update(
                band_lambda_re=float(lb.real),
                band_lambda_f_hz=float(lb.imag / (2 * math.pi)),
                band_pencil_re=float(mb.real),
                band_pencil_f_hz=float(abs(mb.imag) / (2 * math.pi)),
            )
    # a real mode need not move machines against each other: its observable is the
    # deviation of the state with the largest participation, fitted until that
    # deviation reaches 0.05 (its own units)
    left = np.linalg.inv(np.linalg.eig(case.system.A)[1]).T
    values = np.linalg.eigvals(case.system.A)
    c = int(np.argmin(np.abs(values - lam)))
    part = np.abs(np.linalg.eig(case.system.A)[1][:, c] * left[:, c])
    k = int(np.argmax(part))
    s2 = x[k] - case.equilibrium.x[k]
    row["critical_state"] = labels[k]
    aperiodic = abs(lam.imag) / (2 * math.pi) < 1e-3
    if aperiodic:
        small = np.maximum.accumulate(np.abs(s2)) <= 0.05
        t2 = float(times[small][-1]) if small.any() else 0.0
        w2 = small & (times >= t2 / 3.0) & (np.abs(s2) > 0)
        if w2.sum() > 20:
            row["aperiodic_growth_fit"] = float(
                np.polyfit(times[w2], np.log(np.abs(s2[w2])), 1)[0]
            )
            row["aperiodic_window_s"] = [float(times[w2][0]), float(times[w2][-1])]
    row["primary_re"] = (
        row.get("aperiodic_growth_fit") if aperiodic else row.get("pencil_re")
    )
    row["primary_f_hz"] = 0.0 if aperiodic else row.get("pencil_f_hz")
    if diverged:
        outcome = "DIVERGED"
    elif row["env_growth_2nd_half"] > 0:
        outcome = "UNSTABLE"
    else:
        outcome = "STABLE"
    row["outcome"] = outcome
    row["verdict_agrees"] = (outcome in ("DIVERGED", "UNSTABLE")) == bool(rhp)
    stride = max(1, n // 1500)
    trace = pd.DataFrame(
        {"t": times[::stride], "signal": s[::stride], "spread": spread[::stride]}
    )
    return row, trace


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    exp = GateExperiment(
        name="G2_tds",
        question="Does the nonlinear phasor-domain model do what the eigenvalues say?",
        config={
            "kind": "nonlinear phasor-domain DAE, NOT EMT",
            "step": STEP,
            "diverge_pu": DIVERGE,
            "linear_window_pu": LINEAR,
            "horizon": [T_MIN, T_MAX],
            "cases": [c[1] for c in CASES],
            "disturbances": DISTURBANCES,
        },
        workers=8,
    )
    started = time.time()
    tasks = [(*c, d) for c in CASES for d in ("D1", "D2")]
    with Pool(8, initializer=pin_blas_threads) as pool:
        out = pool.map(run, tasks, chunksize=1)
    rows = pd.DataFrame([r for r, _ in out])
    for r, trace in out:
        key = (
            f"{r['case']}_{r['disturbance']}".replace(" ", "_")
            .replace("+", "p")
            .replace("%", "pct")
        )
        trace.to_csv(OUT / f"trace_{key}.csv.gz", index=False)
    rows["f_err_hz"] = (rows.primary_f_hz - rows.lambda_f_hz).abs()
    rows["re_err"] = (rows.primary_re - rows.lambda_re).abs()
    rows["f_err_fft_hz"] = (rows.fft_f_hz - rows.lambda_f_hz).abs()
    rows["growth_sign_agrees"] = np.sign(rows.primary_re) == np.sign(rows.lambda_re)
    rows.to_csv(OUT / "G2_tds_summary.csv", index=False)
    fitted = rows.dropna(subset=["primary_re"])
    summary = {
        "runs": int(len(rows)),
        "runs_with_modal_fit": int(len(fitted)),
        "verdict_agreement": float(rows.verdict_agrees.mean()),
        "growth_sign_agreement": float(fitted.growth_sign_agrees.mean()),
        "median_f_err_hz": float(fitted.f_err_hz.median()),
        "max_f_err_hz": float(fitted.f_err_hz.max()),
        "median_re_err": float(fitted.re_err.median()),
        "max_re_err": float(fitted.re_err.max()),
        "outcomes": rows.groupby(["case", "disturbance"]).outcome.first().to_dict(),
    }
    summary["outcomes"] = {f"{a} | {b}": v for (a, b), v in summary["outcomes"].items()}
    (OUT / "G2_summary.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", elapsed_s=time.time() - started)
    pd.set_option("display.width", 260)
    print(
        rows[
            [
                "case",
                "disturbance",
                "lambda_re",
                "lambda_f_hz",
                "predicted",
                "outcome",
                "primary_re",
                "primary_f_hz",
                "fft_f_hz",
                "env_growth_2nd_half",
                "run_status",
                "t_end_s",
                "cpu_s",
            ]
        ].to_string(index=False)
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "outcomes"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
