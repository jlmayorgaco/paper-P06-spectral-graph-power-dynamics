"""Nonlinear phasor-domain TDS with frozen model-scope guards (final campaign).

NOT EMT. The L0 DAE of a solved case is integrated with the G2 integrator:

- BDF;
- the network solved by chord / Newton at every right-hand side;
- the Jacobian of the reduced field from the DAE blocks at the current point.

A disturbance is a pulse: two segments, and the solver never steps across the
switch. Every run returns exactly one label, following
configs/ias2026/final_nonlinear_composability_v1.yaml:

    RECOVERS / FAILS_DECLARED_SECURITY / OUTSIDE_MODEL_SCOPE / NUMERICAL_FAILURE

The guards stand for limiters and protections that the frozen model omits. The
first violation in time decides OUTSIDE_MODEL_SCOPE, unless a failure was
declared earlier. A trajectory that needs an omitted limiter is never counted as
evidence about composability.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from scipy.integrate import BDF
from scipy.linalg import lu_factor, lu_solve

from ..dynamics.linearize import central_difference_jacobians
from ..models.ieee39_devices import OMEGA_B

F0 = 60.0
LIMITS = Path(__file__).resolve().parents[3] / "configs" / "ias2026" / "ieee39_documented_limits_v1.json"

RECOVERS = "RECOVERS"
FAILS = "FAILS_DECLARED_SECURITY"
OUTSIDE = "OUTSIDE_MODEL_SCOPE"
NUMERICAL = "NUMERICAL_FAILURE"


@dataclass
class Guards:
    bus_v: tuple = (0.80, 1.20)
    gfl_v: tuple = (0.90, 1.10)
    gfl_i: float = 1.20
    gfl_s: float = 1.20
    gfl_q: float = 1.00
    pll_hz: tuple = (-3.0, 1.8)
    speed_hz: tuple = (-3.0, 1.8)
    coi_hz: tuple = (-0.5, 0.5)
    residual: float = 1e-6
    efd: dict = field(default_factory=dict)  # bus -> (vrmin, vrmax)
    pss: dict = field(default_factory=dict)  # bus -> (lsmin, lsmax)

    @classmethod
    def documented(cls, path: Path = LIMITS) -> Guards:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(
            efd={m["bus"]: (m["VRMIN"], m["VRMAX"]) for m in data["machines"]},
            pss={m["bus"]: (m["LSMIN"], m["LSMAX"]) for m in data["machines"]},
        )


@dataclass
class Rules:
    pulse_s: float = 0.5
    after_s: float = 30.0
    max_step: float = 0.05
    window_s: float = 5.0
    growth_factor: float = 2.0
    growth_fit_s: float = 10.0
    decay_factor: float = 0.1
    tail_s: float = 15.0
    s1_scale: float = 1e-3
    s2_scale: float = 1e-2


class _Observer:
    """Observables and guard values of the frozen model at (x, z)."""

    def __init__(self, dae, x_eq, z_eq, guards: Guards):
        self.dae, self.g = dae, guards
        v0 = dae.voltages(z_eq)
        self.v0 = np.abs(v0)
        self.sg, self.gfl = [], []
        for slot in dae.slots:
            p = getattr(slot.device, "parameters", None)
            if slot.kind == "sg" and p is not None and hasattr(p, "m"):
                self.sg.append((slot, p.m * slot.device.weight, p))
            elif slot.kind == "gfl":
                self.gfl.append((slot, slot.device.parameters))
        self.m = np.array([m for _, m, _ in self.sg])
        d_eq = np.array([x_eq[s.start] for s, _, _ in self.sg])
        self.rel0 = d_eq - (self.m @ d_eq) / self.m.sum()

    def __call__(self, x, z) -> tuple[dict, str | None]:
        v = self.dae.voltages(z)
        vm = np.abs(v)
        om = np.array([x[s.start + 1] for s, _, _ in self.sg])
        de = np.array([x[s.start] for s, _, _ in self.sg])
        coi = (self.m @ om) / self.m.sum()
        rel = de - (self.m @ de) / self.m.sum() - self.rel0
        obs = {
            "s1": float(np.abs(om - coi).max()),
            "s2": float(np.abs(vm - self.v0).max()),
            "max_rel_angle": float(np.abs(rel).max()),
            "coi_hz": float((coi - 1.0) * F0),
            "vmin": float(vm.min()),
            "vmax": float(vm.max()),
        }
        viol = None
        if not (self.g.bus_v[0] <= obs["vmin"] and obs["vmax"] <= self.g.bus_v[1]):
            viol = "bus_voltage"
        sp = (om - 1.0) * F0
        obs["speed_hz_min"], obs["speed_hz_max"] = float(sp.min()), float(sp.max())
        if viol is None and not (self.g.speed_hz[0] <= sp.min() and sp.max() <= self.g.speed_hz[1]):
            viol = "machine_speed"
        if viol is None and not (self.g.coi_hz[0] <= obs["coi_hz"] <= self.g.coi_hz[1]):
            viol = "coi_frequency"
        efd_margin, pss_max = np.inf, 0.0
        for slot, _, p in self.sg:
            xs = x[slot.start : slot.stop]
            if slot.bus in self.g.efd:
                lo, hi = self.g.efd[slot.bus]
                efd_margin = min(efd_margin, xs[4] - lo, hi - xs[4])
            out = p.pss_gain * xs[6]
            pss_max = max(pss_max, abs(out))
            if viol is None and slot.bus in self.g.pss:
                lo, hi = self.g.pss[slot.bus]
                if not (lo <= out <= hi):
                    viol = "pss_output"
        obs["efd_margin"], obs["pss_max"] = float(efd_margin), float(pss_max)
        if viol is None and efd_margin < 0:
            viol = "exciter_field_voltage"
        gv, gi, gs, gq, pll = [1.0, 1.0], 0.0, 0.0, 0.0, [0.0, 0.0]
        for slot, p in self.gfl:
            xs = x[slot.start : slot.stop]
            vb = complex(v[self.dae.network.position(slot.bus)])
            rot = vb * np.exp(-1j * xs[0])
            i_mag = math.hypot(xs[6], xs[7])
            q = rot.imag * xs[6] - rot.real * xs[7]
            f_pll = (p.kp_pll * rot.imag + xs[1]) / (2 * math.pi)
            gv = [min(gv[0], abs(vb)), max(gv[1], abs(vb))]
            gi, gs, gq = max(gi, i_mag), max(gs, abs(vb) * i_mag), max(gq, abs(q))
            pll = [min(pll[0], f_pll), max(pll[1], f_pll)]
        obs.update(gfl_vmin=gv[0], gfl_vmax=gv[1], gfl_i=gi, gfl_s=gs, gfl_q=gq,
                   pll_hz_min=pll[0], pll_hz_max=pll[1])
        if viol is None and self.gfl:
            if not (self.g.gfl_v[0] <= gv[0] and gv[1] <= self.g.gfl_v[1]):
                viol = "gfl_terminal_voltage"
            elif gi > self.g.gfl_i:
                viol = "gfl_current"
            elif gs > self.g.gfl_s:
                viol = "gfl_apparent_power"
            elif gq > self.g.gfl_q:
                viol = "gfl_reactive"
            elif not (self.g.pll_hz[0] <= pll[0] and pll[1] <= self.g.pll_hz[1]):
                viol = "pll_frequency"
        return obs, viol


class _Network:
    """z(x) by chord iterations on a frozen LU, then full Newton."""

    def __init__(self, dae, x, z):
        self.z = z.copy()
        self.lu = lu_factor(central_difference_jacobians(dae, x, z, {}).gz)

    def solve(self, dae, x, tol=1e-10):
        z = self.z.copy()
        for _ in range(15):
            r = dae.g(x, z, {})
            if not np.all(np.isfinite(r)):
                break
            if np.abs(r).max() < tol:
                self.z = z
                return z, float(np.abs(r).max())
            z = z - lu_solve(self.lu, r)
        z = self.z.copy()
        for _ in range(12):
            r = dae.g(x, z, {})
            if not np.all(np.isfinite(r)):
                return None, float("inf")
            if np.abs(r).max() < tol:
                self.z = z
                return z, float(np.abs(r).max())
            try:
                self.lu = lu_factor(central_difference_jacobians(dae, x, z, {}).gz)
            except (np.linalg.LinAlgError, ValueError):
                return None, float("inf")
            z = z - lu_solve(self.lu, r)
        r = dae.g(x, z, {})
        return None, float(np.abs(r).max()) if np.all(np.isfinite(r)) else float("inf")


def _envelope_growth(t, y):
    from scipy.signal import find_peaks

    y = np.maximum(np.asarray(y), 1e-300)
    pk, _ = find_peaks(y)
    if pk.size >= 4:
        return float(np.polyfit(t[pk], np.log(y[pk]), 1)[0])
    return float(np.polyfit(t, np.log(y), 1)[0])


def simulate(case, pulse_dae, *, guards: Guards | None = None, rules: Rules | None = None,
             keep_trace: bool = True) -> dict:
    """One pulse run from the equilibrium of ``case``; one label."""

    guards = guards or Guards.documented()
    rules = rules or Rules()
    dae0 = case.dae
    x0, z0 = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    obs_fn = _Observer(dae0, x0, z0, guards)
    net = _Network(dae0, x0, z0)
    trace, label, reason, t_label = [], None, "", None
    peak, post = 0.0, []
    t_pe = rules.pulse_s
    seg_list = [(0.0, rules.pulse_s, pulse_dae), (rules.pulse_s, rules.pulse_s + rules.after_s, dae0)]
    x = x0

    for t0, t1, dae in seg_list:
        def rhs(t, xx, dae=dae):
            z, res = net.solve(dae, xx)
            if z is None:
                return np.full(xx.shape, np.nan)
            return dae.f(xx, z, {})

        def jac(t, xx, dae=dae):
            z, _ = net.solve(dae, xx)
            if z is None:
                return case.system.A
            j = central_difference_jacobians(dae, xx, z, {})
            return j.fx - j.fz @ np.linalg.solve(j.gz, j.gx)

        solver = BDF(rhs, t0, x, t1, jac=jac, rtol=1e-7, atol=1e-9, max_step=rules.max_step)
        while solver.status == "running":
            try:
                solver.step()
            except (np.linalg.LinAlgError, ValueError):
                label, reason, t_label = NUMERICAL, "solver_exception", solver.t
                break
            if solver.status == "failed" or not np.all(np.isfinite(solver.y)):
                label, reason, t_label = NUMERICAL, "solver_failed", solver.t
                break
            z, res = net.solve(dae, solver.y)
            if z is None or res > guards.residual:
                # the network has no nearby solution: voltage collapse or numerics
                label, reason, t_label = NUMERICAL, f"network_residual_{res:.1e}", solver.t
                break
            obs, viol = obs_fn(solver.y, z)
            dval = max(obs["s1"] / rules.s1_scale, obs["s2"] / rules.s2_scale)
            obs.update(t=solver.t, D=dval, residual=res)
            if keep_trace:
                trace.append(obs)
            if viol is not None:
                label, reason, t_label = OUTSIDE, viol, solver.t
                break
            if solver.t >= t_pe:
                post.append((solver.t, dval))
                if solver.t <= t_pe + rules.window_s:
                    peak = max(peak, dval)
                elif peak > 0 and dval > rules.growth_factor * peak:
                    arr = np.array([p for p in post if p[0] >= solver.t - rules.growth_fit_s])
                    if len(arr) > 5 and _envelope_growth(arr[:, 0], arr[:, 1]) > 0:
                        label, reason, t_label = FAILS, "growth", solver.t
                        break
                if obs["max_rel_angle"] > math.pi:
                    label, reason, t_label = FAILS, "loss_of_synchronism", solver.t
                    break
        if label is not None:
            break
        x = solver.y.copy()

    import pandas as pd

    tr = pd.DataFrame(trace)
    if label is None:
        tail = tr[tr.t >= tr.t.iloc[-1] - rules.tail_s]
        final = float(tr.D.iloc[-1])
        growth = _envelope_growth(tail.t.to_numpy(), tail.D.to_numpy()) if len(tail) > 5 else np.nan
        if peak > 0 and final <= rules.decay_factor * peak and growth < 0:
            label, reason = RECOVERS, "decayed"
        else:
            label, reason = FAILS, f"not_recovered final/peak={final / max(peak, 1e-300):.3g} growth={growth:.3g}"
        t_label = float(tr.t.iloc[-1])
    summary = {"label": label, "reason": reason, "t_label": t_label, "D_peak": peak}
    if len(tr):
        for k in ("vmin", "gfl_vmin", "pss_max", "coi_hz", "max_rel_angle", "gfl_i", "efd_margin"):
            if k in tr:
                summary[f"extreme_{k}"] = float(tr[k].min() if k in ("vmin", "gfl_vmin", "efd_margin")
                                                else tr[k].abs().max())
    return {"summary": summary, "trace": tr}


# ------------------------------------------------------------ disturbances --


def pulse_dae(case, family: str, amplitude: float, spec: dict):
    """The DAE during the pulse for one declared disturbance family."""

    from .model import PhasorModel

    dae = case.dae
    if family == "D1":
        model = PhasorModel(case)
        u = np.zeros(model.n_u)
        u[model.input_names.index(f"dP_load_{spec['bus']}")] = amplitude / 100.0
        return model._with_inputs(u)
    if family == "D2":
        model = PhasorModel(case)
        u = np.zeros(model.n_u)
        u[model.input_names.index(f"dPm_sg{spec['machine_bus']}")] = -amplitude / 100.0
        return model._with_inputs(u)
    if family == "D3":
        slots = []
        for slot in dae.slots:
            if slot.kind == "gfl":
                p = slot.device.parameters
                dev = replace(slot.device, parameters=replace(p, p_ref=p.p_ref * (1.0 - amplitude)))
                slot = replace(slot, device=dev)
            slots.append(slot)
        return replace(dae, slots=tuple(slots))
    raise ValueError(family)


def transverse_alpha(case) -> float:
    from ..certification.symmetry import frequency_partner, rotation_generator
    from ..certification.transverse import transverse_operator

    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    return float(np.linalg.eigvals(transverse_operator(case.system.A, r_x, w).a_perp).real.max())


__all__ = ["simulate", "pulse_dae", "Guards", "Rules", "transverse_alpha"]
