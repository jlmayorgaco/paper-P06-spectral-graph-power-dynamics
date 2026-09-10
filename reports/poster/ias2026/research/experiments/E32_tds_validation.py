"""E32 - overnight 3. Nonlinear phasor-domain DAE time-domain validation.

This is NOT an electromagnetic-transient study. It is a nonlinear
phasor-domain DAE simulation: the network is represented by its algebraic
power-flow equations and the devices by their differential equations, integrated
without linearization. It tests whether the modes the eigenanalysis reports are
the modes the nonlinear system actually exhibits.

Six configurations, four disturbances, identical magnitude everywhere:

    D1  +2 % active load at bus 20
    D2  +2 % active load at bus 29
    D3  +2 % mechanical power at the machine on bus 31, present in every case
    D4  +2 % active-power reference at the converter on bus 30, converter cases only

The step is applied as an initial-condition mismatch rather than as a
discontinuity inside the right-hand side, which is the formulation E25 arrived at
after the stiff solver stalled on the switched version.

Frequency is estimated two independent ways, an FFT peak and a matrix-pencil
ringdown fit, and both are compared against the eigenvalue prediction.
"""

from __future__ import annotations

import time
from dataclasses import replace

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import lu_factor, lu_solve

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from E25_time_domain_validation import match_mode, matrix_pencil
from ibr_cycles.dynamics.linearize import central_difference_jacobians
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

NAME = "E32_TDS_validation"
STEP = 0.02
T_END = 12.0
T_SETTLE = 0.3
LINEARITY_BOUND = 2.0e-4
PM_BUS = 31          # a machine present in every configuration
#: A 2 % mechanical step on one machine is not small: it drives the speed spread
#: past the small-signal window before a ringdown can be fitted. D3 therefore
#: carries its own magnitude, identical across every configuration, which is what
#: the protocol requires. It is declared here, not tuned per case.
PM_STEP = 0.002
PREF_BUS = 30        # a converter, when the configuration has one

#: RC, reconstructed from the frozen E21 record. tau_p was not among the
#: reported multipliers, so it is taken unchanged and the reconstruction is
#: verified against the recorded abscissa of -0.0500 before it is used.
RC_MULTIPLIERS = {
    "kp_p": 1.226,
    "ki_p": 0.460,
    "kp_q": 1.173,
    "ki_q": 1.273,
    "kp_i": 1.115,
    "ki_i": 0.960,
}
RB_SCALE = 0.4273112886362242


def rc_converter():
    base = ConverterParameters()
    return replace(base, **{k: getattr(base, k) * v for k, v in RC_MULTIPLIERS.items()})


def rb_converter():
    base = ConverterParameters()
    return replace(base, kp_p=base.kp_p * RB_SCALE, ki_p=base.ki_p * RB_SCALE)


def perturbed_devices(case, kind):
    """Return a slot tuple with one device setpoint stepped, or None."""

    slots = list(case.dae.slots)
    for position, slot in enumerate(slots):
        if kind == "pm" and slot.kind == "sg" and slot.bus == PM_BUS:
            device = slot.device
            parameters = replace(
                device.parameters, pm=device.parameters.pm * (1.0 + PM_STEP)
            )
            slots[position] = replace(
                slot, device=replace(device, parameters=parameters)
            )
            return tuple(slots)
        if kind == "pref" and slot.kind == "gfl" and slot.bus == PREF_BUS:
            device = slot.device
            parameters = replace(
                device.parameters, p_ref=device.parameters.p_ref * (1.0 + STEP)
            )
            slots[position] = replace(
                slot, device=replace(device, parameters=parameters)
            )
            return tuple(slots)
    return None


def simulate(case, disturbance):
    """Integrate the nonlinear DAE from a stepped setpoint or load."""

    dae = case.dae
    network = dae.network
    x0, z0 = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    factor = lu_factor(central_difference_jacobians(dae, x0, z0, {}).gz)
    reduced = case.system.A
    speeds = [i for i, n in enumerate(case.system.labels) if n.startswith("omega_sg")]

    original_load = None
    original_slots = None
    if disturbance["kind"] == "load":
        bus = disturbance["bus"]
        original_load = network.loads[bus]
        network.loads[bus] = complex(
            original_load.real * (1.0 + STEP), original_load.imag
        )
    else:
        stepped = perturbed_devices(case, disturbance["kind"])
        if stepped is None:
            return None, 0
        original_slots = dae.slots
        object.__setattr__(dae, "slots", stepped)

    state = {"z": z0.copy(), "factor": factor}

    def rhs(t, x):
        z = state["z"]
        for attempt in range(2):
            for _ in range(8):
                residual = dae.g(x, z, {})
                if np.abs(residual).max() < 1e-11:
                    break
                z = z - lu_solve(state["factor"], residual)
            else:
                if attempt == 0:
                    state["factor"] = lu_factor(
                        central_difference_jacobians(dae, x, z, {}).gz
                    )
                    continue
            break
        state["z"] = z
        return dae.f(x, z, {})

    def leaves_small_signal(t, x):
        spread = x[speeds] - x[speeds].mean()
        return LINEARITY_BOUND - float(np.abs(spread).max())

    leaves_small_signal.terminal = True
    leaves_small_signal.direction = -1

    try:
        solution = solve_ivp(
            rhs,
            (0.0, T_END),
            x0,
            method="BDF",
            jac=lambda t, x: reduced,
            dense_output=True,
            rtol=1e-7,
            atol=1e-9,
            events=leaves_small_signal,
        )
    finally:
        if original_load is not None:
            network.loads[disturbance["bus"]] = original_load
        if original_slots is not None:
            object.__setattr__(dae, "slots", original_slots)
    return solution, len(solution.t)


def inter_area_signal(case, solution):
    labels = case.system.labels
    speeds = [i for i, n in enumerate(labels) if n.startswith("omega_sg")]
    if len(speeds) < 2:
        return None, None
    stop = float(solution.t[-1])
    if stop <= T_SETTLE + 0.6:
        return None, None
    times = np.linspace(T_SETTLE, stop, 4000)
    values = solution.sol(times)
    return times, values[speeds[0]] - values[speeds[-1]]


def fft_peak(times, signal, band=BAND_HZ):
    y = np.asarray(signal) - np.mean(signal)
    window = np.hanning(y.size)
    spectrum = np.abs(np.fft.rfft(y * window))
    frequencies = np.fft.rfftfreq(y.size, d=float(times[1] - times[0]))
    keep = (frequencies >= band[0]) & (frequencies <= band[1])
    if not keep.any():
        return float("nan")
    return float(frequencies[keep][int(np.argmax(spectrum[keep]))])


def envelope_growth(times, signal):
    """Least-squares growth rate of the absolute-value envelope, sign only."""

    y = np.abs(np.asarray(signal))
    y = np.maximum(y, 1e-18)
    slope = np.polyfit(times, np.log(y), 1)[0]
    return float(slope)


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Do the nonlinear trajectories carry the modes the eigenanalysis predicts?",
        config={
            "kind": "nonlinear phasor-domain DAE TDS, NOT EMT",
            "step_fraction": STEP,
            "t_end": T_END,
            "linearity_bound": LINEARITY_BOUND,
            "disturbances": {
                "D1": "+2% active load at bus 20",
                "D2": "+2% active load at bus 29",
                "D3": f"+{PM_STEP:.1%} mechanical power at the machine on bus 31",
                "D4": "+2% active-power reference at the converter on bus 30",
            },
            "rc_multipliers": RC_MULTIPLIERS,
            "rb_scale": RB_SCALE,
        },
        workers=1,
    )
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship_plan = ReplacementPlan.of({b: 1.0 for b in CORE})
    flagship = solve_case(flagship_plan)
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)

    rc = rc_converter()
    rc_case = solve_case(flagship_plan, converter=rc)
    rc_abscissa = max(
        m.real for m in eigen_analysis(rc_case.system.A).modes if abs(m.value) > 1e-3
    )
    experiment.note(
        f"RC reconstruction check: spectral abscissa {rc_abscissa:+.4f} "
        f"against the frozen E21 record of -0.0500"
    )
    rc_ok = abs(rc_abscissa + 0.05) < 0.02

    configurations = {
        "base": base,
        "closest stable proper subset (37)": solve_case(
            ReplacementPlan.of({37: 1.0})
        ),
        "flagship 30+33+35+37": flagship,
        "RC minimum-norm retune": rc_case,
        "RB weaken P loop": solve_case(flagship_plan, converter=rb_converter()),
        "SC condenser 25%": solve_case(
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            )
        ),
        "SG restoration (keep 30)": solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE if b != 30})
        ),
    }
    disturbances = [
        {"name": "D1 load bus 20", "kind": "load", "bus": 20},
        {"name": "D2 load bus 29", "kind": "load", "bus": 29},
        {"name": "D3 pm bus 31", "kind": "pm"},
        {"name": "D4 pref bus 30", "kind": "pref"},
    ]

    rows = []
    for label, case in configurations.items():
        spectrum = eigen_analysis(case.system.A)
        reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
        predicted = reading.family.worst.value if reading else complex("nan")
        alpha = reading.family.alpha if reading else float("nan")
        frequency = reading.family.frequency_worst_hz if reading else float("nan")
        for disturbance in disturbances:
            row = {
                "configuration": label,
                "disturbance": disturbance["name"],
                "alpha_predicted": alpha,
                "freq_predicted_hz": frequency,
            }
            try:
                solution, steps = simulate(case, disturbance)
            except Exception as error:  # noqa: BLE001
                row["status"] = f"ERROR: {str(error)[:80]}"
                rows.append(row)
                continue
            if solution is None:
                row["status"] = "NOT APPLICABLE"
                rows.append(row)
                continue
            times, signal = inter_area_signal(case, solution)
            row["t_end_reached"] = float(solution.t[-1])
            row["terminated_early"] = bool(solution.status == 1)
            if times is None:
                row["status"] = "WINDOW TOO SHORT"
                rows.append(row)
                continue
            exponents = matrix_pencil(times, signal, order=14)
            matched, distance = match_mode(exponents, predicted)
            row.update(
                {
                    "status": "OK",
                    "fft_freq_hz": fft_peak(times, signal),
                    "pencil_freq_hz": float(abs(matched.imag) / (2 * np.pi)),
                    "pencil_alpha": float(matched.real),
                    "pencil_distance": distance,
                    "envelope_growth": envelope_growth(times, signal),
                    "signal_peak": float(np.abs(signal).max()),
                }
            )
            row["freq_error_fft_hz"] = abs(row["fft_freq_hz"] - frequency)
            row["freq_error_pencil_hz"] = abs(row["pencil_freq_hz"] - frequency)
            row["growth_sign_agrees"] = bool(
                np.sign(row["envelope_growth"]) == np.sign(alpha)
            )
            rows.append(row)
            print(
                "  %-34s %-16s pred %+.4f @ %.4f Hz | fft %.4f pencil %.4f a=%+.4f"
                % (
                    label,
                    disturbance["name"],
                    alpha,
                    frequency,
                    row["fft_freq_hz"],
                    row["pencil_freq_hz"],
                    row["pencil_alpha"],
                ),
                flush=True,
            )

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E32_TDS_summary.csv", parquet=True)
    good = table[table.status == "OK"]
    validation = good[
        [
            "configuration",
            "disturbance",
            "alpha_predicted",
            "freq_predicted_hz",
            "fft_freq_hz",
            "pencil_freq_hz",
            "pencil_alpha",
            "freq_error_fft_hz",
            "freq_error_pencil_hz",
            "growth_sign_agrees",
        ]
    ]
    experiment.save_table(validation, "E32_TDS_frequency_validation.csv")

    checks = {
        "runs_ok": int(len(good)),
        "runs_total": int(len(table)),
        "rc_reconstruction_ok": bool(rc_ok),
        "median_fft_freq_error_hz": float(good.freq_error_fft_hz.median()),
        "max_fft_freq_error_hz": float(good.freq_error_fft_hz.max()),
        "median_pencil_freq_error_hz": float(good.freq_error_pencil_hz.median()),
        "growth_sign_agreement": float(good.growth_sign_agrees.mean()),
        "unstable_cases_grow": bool(
            good[good.alpha_predicted > 0].envelope_growth.gt(0).all()
        )
        if (good.alpha_predicted > 0).any()
        else None,
        "stable_cases_decay": bool(
            good[good.alpha_predicted < 0].envelope_growth.lt(0).all()
        ),
    }
    verdict = (
        "PASS"
        if checks["median_fft_freq_error_hz"] < 0.05
        and checks["growth_sign_agreement"] >= 0.9
        else "MIXED"
    )
    print()
    for key, value in checks.items():
        print(f"  {key:32s} {value}")
    experiment.finish(verdict, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
