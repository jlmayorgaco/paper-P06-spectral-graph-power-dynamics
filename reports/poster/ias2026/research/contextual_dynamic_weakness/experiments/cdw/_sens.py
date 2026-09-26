# ruff: noqa: E501
"""CDW sensitivity engine: equilibrium semantics (SPR / RP), IFT total derivatives,
frozen partials and finite re-equilibrated validation (docs/CDW_THEORY_V1.md C4).

w = (x, z, r_free). R(w; a) = [f; g; schedule equations]. Everything uses the TX4
device and network code; only references and parameters are swapped in.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import scipy.linalg as sla

import _cdw as C
from ibr_cycles.dynamics.linearize import central_difference_jacobians
from ibr_cycles.models.ieee39_case import DeviceSlot, Ieee39Dae

SG_REFS = ("pm", "vref")
GFL_REFS = ("p_ref", "q_ref", "v_ref")


def _kind(slot):
    name = type(slot.device).__name__
    return {"SynchronousMachine": "sg", "GridFollowingConverter": "gfl", "StaticInjection": "static"}[name]


class Engine:
    """One solved portfolio at one policy (and optional envelope draw)."""

    def __init__(self, members, theta, *, draw=None, device="gfl", design=None):
        self.members = tuple(members)
        self.theta = tuple(theta)
        self.draw = draw
        self.device = device
        self.design = design or {}
        self.line_scale = {int(e): float(v) for e, v in (self.design.get("line") or {}).items()}
        self.case = C.solve(self.members, self.theta, draw=draw, device=device, **design_kwargs(self.design))
        dae = self.case.dae
        self.net0 = dae.network
        self.pf = dae.power_flow
        self.plan = dae.plan
        self.slots0 = list(dae.slots)
        self.kinds = [_kind(s) for s in self.slots0]
        self.x0 = self.case.equilibrium.x.copy()
        self.z0 = self.case.equilibrium.z.copy()
        self.nx, self.nz = self.x0.size, self.z0.size
        self.slack = self.net0.slack_bus
        self.pos = {s.bus: self.net0.position(s.bus) for s in self.slots0}
        self.r0 = self._refs_of(self.slots0)

    # ------------------------------------------------------------ references --
    def _refs_of(self, slots):
        r = []
        for s, k in zip(slots, self.kinds, strict=True):
            p = s.device.parameters if k != "static" else None
            if k == "sg":
                r += [p.pm, p.vref]
            elif k == "gfl":
                r += [p.p_ref, p.q_ref, p.v_ref]
            else:
                r += [s.device.injection_pu.real, s.device.injection_pu.imag]
        return np.array(r, float)

    def _ref_slices(self):
        out, c = [], 0
        for k in self.kinds:
            n = {"sg": 2, "gfl": 3, "static": 2}[k]
            out.append((c, c + n))
            c += n
        return out

    def free_mask(self, semantics):
        mask = np.zeros(self.r0.size, bool)
        for (a, b), s, k in zip(self._ref_slices(), self.slots0, self.kinds, strict=True):
            if semantics == "SPR":
                mask[a:b] = True
            elif semantics == "RP" and k == "sg" and s.bus == self.slack:
                mask[a] = True  # slack pm only
        return mask

    # ------------------------------------------------------------ parameters --
    def model(self, a: dict | None, r: np.ndarray):
        """Return (dae, schedule) for parameter spec a = {'kind', 'idx', 'value'} and refs r."""

        net = self.net0
        sched_v = {b: self.net0.pv[b]["v"] for b in self.net0.pv}
        sched_v[self.slack] = self.net0.slack_voltage
        sched_p = {b: self.net0.pv[b]["p"] for b in self.net0.pv}
        dev_mod = {}
        rr = r.copy()
        if a:
            kind, idx, val = a["kind"], a["idx"], a["value"]
            if kind == "line":
                net = replace(net, ybus=C.build_ybus({**self.line_scale, idx: self.line_scale.get(idx, 1.0) * val}))
            elif kind == "load":
                net = replace(net, loads={b: v * (val if b == idx else 1.0) for b, v in net.loads.items()})
            elif kind == "vset":
                if a.get("semantics", "SPR") == "SPR":
                    sched_v[idx] = sched_v[idx] + val
                else:  # RP: device voltage reference
                    for (lo, hi), s, k in zip(self._ref_slices(), self.slots0, self.kinds, strict=True):
                        if s.bus == idx:
                            if k == "sg":
                                rr[lo + 1] += val
                            elif k == "gfl":
                                rr[lo + 2] += val
            elif kind in ("g", "pll", "ka"):
                dev_mod[idx] = (kind, val)
        slots = []
        for (lo, hi), s, k in zip(self._ref_slices(), self.slots0, self.kinds, strict=True):
            dev = s.device
            if k == "sg":
                p = replace(dev.parameters, pm=float(rr[lo]), vref=float(rr[lo + 1]))
                if s.bus in dev_mod and dev_mod[s.bus][0] == "ka":
                    p = replace(p, ka=dev.parameters.ka * dev_mod[s.bus][1])
                dev = replace(dev, parameters=p)
            elif k == "gfl":
                p = replace(dev.parameters, p_ref=float(rr[lo]), q_ref=float(rr[lo + 1]), v_ref=float(rr[lo + 2]))
                if s.bus in dev_mod:
                    kk, val = dev_mod[s.bus]
                    if kk == "g":
                        p = replace(p, voltage_gain=float(val))
                    elif kk == "pll":
                        p = replace(p, kp_pll=dev.parameters.kp_pll * val, ki_pll=dev.parameters.ki_pll * val)
                dev = replace(dev, parameters=p)
            else:
                dev = replace(dev, injection_pu=complex(rr[lo], rr[lo + 1]))
            slots.append(DeviceSlot(device=dev, start=s.start, stop=s.stop, kind=s.kind, bus=s.bus, weight=s.weight, loading=s.loading))
        dae = Ieee39Dae(network=net, power_flow=self.pf, slots=tuple(slots), plan=self.plan)
        return dae, (sched_p, sched_v)

    # --------------------------------------------------------------- residual --
    def residual(self, w, a, semantics):
        mask = self.free_mask(semantics)
        x, z = w[: self.nx], w[self.nx: self.nx + self.nz]
        r = self.r0.copy()
        r[mask] = w[self.nx + self.nz:]
        dae, (sp, sv) = self.model(a, r)
        res = [dae.f(x, z, {}), dae.g(x, z, {})]
        v = z[0::2] + 1j * z[1::2]
        extra = []
        for (lo, hi), s, k in zip(self._ref_slices(), dae.slots, self.kinds, strict=True):
            vi = v[self.pos[s.bus]]
            if semantics == "SPR":
                s_dev = vi * np.conj(s.device.injection(x[s.start: s.stop], complex(vi)))
                if k == "sg" and s.bus == self.slack:
                    extra += [abs(vi) - sv[s.bus], np.angle(vi) - self.net0.slack_angle]
                elif k == "sg" or k == "static":
                    extra += [s_dev.real - sp[s.bus], abs(vi) - sv[s.bus]]
                else:  # gfl
                    extra += [s_dev.real - sp[s.bus], abs(vi) - sv[s.bus], r[lo + 2] - abs(vi)]
            elif semantics == "RP" and k == "sg" and s.bus == self.slack:
                extra += [np.angle(vi) - self.net0.slack_angle]
        return np.concatenate(res + [np.array(extra, float)])

    def w0(self, semantics):
        return np.concatenate([self.x0, self.z0, self.r0[self.free_mask(semantics)]])

    def jac_R(self, w, a, semantics, rel=1e-7):
        n = w.size
        J = np.zeros((n, n))
        for j in range(n):
            h = rel * max(1.0, abs(w[j]))
            wp, wm = w.copy(), w.copy()
            wp[j] += h
            wm[j] -= h
            J[:, j] = (self.residual(wp, a, semantics) - self.residual(wm, a, semantics)) / (2 * h)
        return J

    def solve(self, a, semantics, w_init=None, tol=1e-10, maxit=25):
        w = self.w0(semantics) if w_init is None else w_init.copy()
        for it in range(maxit):
            R = self.residual(w, a, semantics)
            nr = float(np.abs(R).max())
            if nr < tol:
                return w, nr, it
            J = self.jac_R(w, a, semantics)
            step = np.linalg.solve(J, -R)
            lam = 1.0
            for _ in range(8):
                wn = w + lam * step
                if float(np.abs(self.residual(wn, a, semantics)).max()) < nr:
                    break
                lam *= 0.5
            w = wn
        R = self.residual(w, a, semantics)
        return w, float(np.abs(R).max()), maxit

    # ------------------------------------------------------------ linear model --
    def A_of(self, w, a, semantics):
        mask = self.free_mask(semantics)
        x, z = w[: self.nx], w[self.nx: self.nx + self.nz]
        r = self.r0.copy()
        r[mask] = w[self.nx + self.nz:]
        dae, _ = self.model(a, r)
        jac = central_difference_jacobians(dae, x, z, {})
        self._last = (dae, z)
        return jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx), jac

    def transverse(self, A, dae=None, z=None):
        """Transverse eigenvalues (TX4 operator) of A at the model/point last built."""

        dae, z = (dae, z) if dae is not None else self._last
        r_x, _ = C.rotation_generator(dae, z)
        w = C.frequency_partner(dae).w
        tr = C.transverse_operator(A, r_x, w)
        return np.linalg.eigvals(tr.a_perp)

    def critical_T(self, A, dae=None, z=None):
        """Rightmost TRANSVERSE eigenvalue, with left/right eigenvectors of the full A."""

        ev = self.transverse(A, dae, z)
        lam_t = complex(ev[np.argmax(ev.real)])
        if lam_t.imag < 0:
            lam_t = lam_t.conjugate()
        vals, vl, vr = sla.eig(A, left=True, right=True)
        j = int(np.argmin(np.abs(vals - lam_t)))
        others = ev[(np.abs(ev - lam_t) > 1e-8) & (np.abs(ev - lam_t.conjugate()) > 1e-8)]
        gap = float(lam_t.real - others.real.max()) if others.size else float("inf")
        return complex(vals[j]), vl[:, j], vr[:, j], gap, float(abs(vals[j] - lam_t))

    @staticmethod
    def critical(A):
        vals, vl, vr = sla.eig(A, left=True, right=True)
        ok = np.abs(vals) > C.STRUCT_ZERO
        idx = np.where(ok)[0]
        j = idx[np.argmax(vals[idx].real)]
        lam = vals[j]
        if lam.imag < 0:  # take the upper-half-plane member
            cands = idx[np.abs(vals[idx] - lam.conjugate()) < 1e-8]
            if cands.size:
                j = cands[0]
                lam = vals[j]
        others = vals[idx]
        others = others[(np.abs(others - lam) > 1e-8) & (np.abs(others - lam.conjugate()) > 1e-8)]
        gap = float(lam.real - others.real.max()) if others.size else float("inf")
        return complex(lam), vl[:, j], vr[:, j], gap

    # ----------------------------------------------------------- derivatives --
    def derivatives(self, specs, semantics, *, h_rel=1e-3):
        """Frozen and total d lambda_c / d a for every spec (base value in spec['value']).

        dw*/da = -R_w^{-1} R_a (IFT). dA_frozen = [A(w0, a+h) - A(w0, a-h)] / 2h;
        dA_total = [A(w0 + h dw, a+h) - A(w0 - h dw, a-h)] / 2h (directional derivative of
        A(w*(a), a)); d lambda = v^H dA u / v^H u.
        """

        w0 = self.w0(semantics)
        R0 = self.residual(w0, None, semantics)
        A0, _ = self.A_of(w0, None, semantics)
        lam, vl, vr, gap, match_err = self.critical_T(A0)
        denom = np.vdot(vl, vr)
        Rw = self.jac_R(w0, None, semantics)
        lu = sla.lu_factor(Rw)
        out = []
        for sp in specs:
            base = sp["value"]
            h = h_rel * max(1.0, abs(base))
            ap = {**sp, "semantics": semantics, "value": base + h}
            am = {**sp, "semantics": semantics, "value": base - h}
            Ra = (self.residual(w0, ap, semantics) - self.residual(w0, am, semantics)) / (2 * h)
            dw = -sla.lu_solve(lu, Ra)
            App, _ = self.A_of(w0, ap, semantics)
            Amm, _ = self.A_of(w0, am, semantics)
            dA_part = (App - Amm) / (2 * h)
            Atp, _ = self.A_of(w0 + h * dw, ap, semantics)
            Atm, _ = self.A_of(w0 - h * dw, am, semantics)
            dA_tot = (Atp - Atm) / (2 * h)
            d_frozen = complex(np.vdot(vl, dA_part @ vr) / denom)
            d_total = complex(np.vdot(vl, dA_tot @ vr) / denom)
            out.append({**{k: v for k, v in sp.items() if k != "semantics"}, "semantics": semantics,
                        "d_frozen": [d_frozen.real, d_frozen.imag], "d_total": [d_total.real, d_total.imag],
                        "dw_norm": float(np.abs(dw).max()), "h": h})
        return {"lam": [lam.real, lam.imag], "gap2": gap, "R0": float(np.abs(R0).max()), "match_err": match_err, "items": out}

    # --------------------------------------------------------------- finite --
    def finite(self, sp, value, semantics, lam_ref):
        """Finite re-equilibrated alpha and tracked eigenvalue at parameter value."""

        if semantics == "SPR":
            case = spr_case(self, sp, value)
            if case is None:
                return None
            A = C.physical_matrices(case)[0]
            ev = transverse_spectrum(case, A)
        else:
            w, nr, _ = self.solve({**sp, "semantics": "RP", "value": value}, "RP")
            if nr > 1e-8:
                return None
            A, _ = self.A_of(w, {**sp, "semantics": "RP", "value": value}, "RP")
            ev = self.transverse(A)
        j = int(np.argmin(np.abs(ev - lam_ref)))
        return {"alpha": float(ev.real.max()), "tracked_re": float(ev[j].real), "tracked_hz": float(abs(ev[j].imag) / (2 * np.pi))}


def transverse_spectrum(case, A):
    r_x, _ = C.rotation_generator(case.dae, case.equilibrium.z)
    w = C.frequency_partner(case.dae).w
    tr = C.transverse_operator(A, r_x, w)
    return np.linalg.eigvals(tr.a_perp)


def design_kwargs(design, extra=None):
    """solve kwargs for a design point {'g': {bus: abs}, 'ka': {bus: factor}, 'line': {e: scale}}
    plus one extra spec change (kind, idx, value) applied on top."""

    g = {int(b): float(v) for b, v in (design.get("g") or {}).items()}
    ka = {int(b): float(v) for b, v in (design.get("ka") or {}).items()}
    line = {int(e): float(v) for e, v in (design.get("line") or {}).items()}
    loads, vset, pll = None, None, {}
    if extra:
        kind, idx, value = extra
        if kind == "line":
            line = {**line, idx: line.get(idx, 1.0) * value}
        elif kind == "load":
            loads = {idx: value}
        elif kind == "vset":
            vset = {idx: value}
        elif kind == "g":
            g = {**g, idx: value}
        elif kind == "ka":
            ka = {**ka, idx: ka.get(idx, 1.0) * value}
        elif kind == "pll":
            pll = {idx: value}
    kw = {}
    conv = {b: {"voltage_gain": ("abs", v)} for b, v in g.items()}
    for b, v in pll.items():
        conv.setdefault(b, {}).update({"kp_pll": v, "ki_pll": v})
    if conv:
        kw["conv_extra"] = conv
    if ka:
        kw["mbs_extra"] = {b: {"ka": v} for b, v in ka.items()}
    if line or loads or vset:
        kw["network"] = C.network_with(scale=line or None, loads_scale=loads, vset=vset)
    return kw


def spr_case(eng: Engine, sp, value):
    """Independent SPR truth: TX4 solve_case with the parameter applied (on top of the design)."""

    kind, idx = sp["kind"], sp["idx"]
    kw = {"draw": eng.draw, "device": eng.device, **design_kwargs(eng.design, (kind, idx, value))}
    try:
        return C.solve(eng.members, eng.theta, **kw)
    except (C.InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
        return None
