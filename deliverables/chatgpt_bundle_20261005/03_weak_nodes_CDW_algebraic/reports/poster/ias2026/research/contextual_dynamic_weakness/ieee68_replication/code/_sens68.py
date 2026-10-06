# ruff: noqa: E501
"""CDW68 sensitivity engine (Model A): SPR-68 semantics and implicit-function total derivatives of the critical
transverse eigenvalue with respect to whole-two-port branch scaling (docs/CDW68_PREREG_V1.md section 11).

w = (x, z, r). R(w; gamma) = [f; g; schedules]:
  surviving machine (non-slack): P_gen = P_sched, |V| = V_sched      freed: P-reference (governor pref or tm), V-reference (vref or manual Efd)
  slack machine G16:             |V| = V_sched, angle(V) = 0         freed: P-reference, manual Efd
  converter:                     P = P_sched, |V| = V_sched, v_ref = |V|   freed: p_ref, q_ref, v_ref
  load:                          S_load(V) = S_sched                  freed: load admittance (G, B)
Everything else is the frozen 68-bus DAE (code/_r68.build68). Port of experiments/cdw/_sens.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import scipy.linalg as sla

import _r68 as R
from ibr_cycles.dynamics.linearize import central_difference_jacobians
from ibr_cycles.models.ieee39_case import Ieee39Dae


@dataclass
class Dae68(Ieee39Dae):
    """The frozen DAE with explicit constant load admittances (same equations as the 68-bus impedance loads)."""

    yload: dict = field(default_factory=dict)

    def g(self, x, z, theta):
        v = self.voltages(z)
        injection = np.zeros(self.network.n_bus, dtype=np.complex128)
        for slot in self.slots:
            pos = self.network.position(slot.bus)
            injection[pos] += slot.device.injection(x[slot.start: slot.stop], complex(v[pos]))
        for bus, yl in self.yload.items():
            pos = self.network.position(bus)
            injection[pos] -= yl * v[pos]
        residual = self.network.ybus @ v - injection
        out = np.empty(self.n_z)
        out[0::2] = residual.real
        out[1::2] = residual.imag
        return out


def _is_gov(dev):
    return isinstance(dev, R.Gov68)


def _mach(dev):
    return dev.base if _is_gov(dev) else dev


class Engine68:
    def __init__(self, members, *, g=0.0, k=1.0, variant="REAL", draw=None):
        self.members = tuple(members)
        self.variant = variant
        self.case = R.build68(self.members, variant=variant, g=g, k=k, draw=draw)
        dae = self.case.dae
        self.net0 = dae.network
        self.pf = dae.power_flow
        self.slots0 = list(dae.slots)
        self.x0 = self.case.equilibrium.x.copy()
        self.z0 = self.case.equilibrium.z.copy()
        self.nx, self.nz = self.x0.size, self.z0.size
        self.slack = self.net0.slack_bus
        self.loads = sorted(self.net0.loads)
        v0 = self.pf.voltages
        self.y0 = {b: np.conj(self.net0.loads[b]) / abs(v0[self.net0.position(b)]) ** 2 for b in self.loads}
        self.pos = {s.bus: self.net0.position(s.bus) for s in self.slots0}
        self.r0 = self._refs_of(self.slots0)
        self._check_dae()

    # ------------------------------------------------------------ references --
    def _refs_of(self, slots):
        r = []
        for s in slots:
            dev = s.device
            if s.kind == "sg":
                p = _mach(dev).parameters
                r += [dev.pref if _is_gov(dev) else p.tm, p.vref if p.exciter in ("DC4B", "ST1A") else p.efd0]
            else:
                p = dev.parameters
                r += [p.p_ref, p.q_ref, p.v_ref]
        for b in self.loads:
            r += [self.y0[b].real, self.y0[b].imag]
        return np.array(r, float)

    def _slices(self):
        out, c = [], 0
        for s in self.slots0:
            n = 2 if s.kind == "sg" else 3
            out.append((c, c + n))
            c += n
        return out, c

    def model(self, a, r):
        net = self.net0
        if a:
            net = R.network_with(scale={int(a["idx"]): float(a["value"])})
        sl, c0 = self._slices()
        slots = []
        for (lo, hi), s in zip(sl, self.slots0, strict=True):
            dev = s.device
            if s.kind == "sg":
                m = _mach(dev)
                p = m.parameters
                pv = float(r[lo + 1])
                p = replace(p, vref=pv) if p.exciter in ("DC4B", "ST1A") else replace(p, efd0=pv)
                if _is_gov(dev):
                    dev = replace(dev, base=replace(m, parameters=p), pref=float(r[lo]))
                else:
                    dev = replace(m, parameters=replace(p, tm=float(r[lo])))
            else:
                p = replace(dev.parameters, p_ref=float(r[lo]), q_ref=float(r[lo + 1]), v_ref=float(r[lo + 2]))
                dev = replace(dev, parameters=p)
            slots.append(replace(s, device=dev))
        yl = {b: complex(r[c0 + 2 * j], r[c0 + 2 * j + 1]) for j, b in enumerate(self.loads)}
        dae = Dae68(network=net, power_flow=self.pf, slots=tuple(slots), plan=self.case.dae.plan, yload=yl)
        return dae

    def _check_dae(self):
        dae = self.model(None, self.r0)
        e = max(np.abs(dae.f(self.x0, self.z0, {})).max(), np.abs(dae.g(self.x0, self.z0, {})).max())
        assert e < 1e-7, f"Dae68 does not reproduce the frozen equilibrium ({e:.2e})"

    # ---------------------------------------------------------------- residual --
    def residual(self, w, a):
        x, z = w[: self.nx], w[self.nx: self.nx + self.nz]
        r = w[self.nx + self.nz:]
        dae = self.model(a, r)
        v = z[0::2] + 1j * z[1::2]
        sl, c0 = self._slices()
        extra = []
        for (lo, hi), s in zip(sl, dae.slots, strict=True):
            vi = v[self.pos[s.bus]]
            sdev = vi * np.conj(s.device.injection(x[s.start: s.stop], complex(vi)))
            if s.kind == "sg" and s.bus == self.slack:
                extra += [abs(vi) - self.net0.slack_voltage, np.angle(vi) - self.net0.slack_angle]
            elif s.kind == "sg":
                extra += [sdev.real - self.net0.pv[s.bus]["p"], abs(vi) - self.net0.pv[s.bus]["v"]]
            else:
                extra += [sdev.real - self.net0.pv[s.bus]["p"], abs(vi) - self.net0.pv[s.bus]["v"], r[lo + 2] - abs(vi)]
        for j, b in enumerate(self.loads):
            vb = v[self.net0.position(b)]
            yl = complex(r[c0 + 2 * j], r[c0 + 2 * j + 1])
            s_load = np.conj(yl) * abs(vb) ** 2
            extra += [s_load.real - self.net0.loads[b].real, s_load.imag - self.net0.loads[b].imag]
        return np.concatenate([dae.f(x, z, {}), dae.g(x, z, {}), np.array(extra, float)])

    def w0(self):
        return np.concatenate([self.x0, self.z0, self.r0])

    def jac_R(self, w, a, rel=1e-7):
        n = w.size
        J = np.zeros((n, n))
        for j in range(n):
            h = rel * max(1.0, abs(w[j]))
            wp, wm = w.copy(), w.copy()
            wp[j] += h
            wm[j] -= h
            J[:, j] = (self.residual(wp, a) - self.residual(wm, a)) / (2 * h)
        return J

    # ----------------------------------------------------------- linear model --
    def A_of(self, w, a):
        x, z = w[: self.nx], w[self.nx: self.nx + self.nz]
        dae = self.model(a, w[self.nx + self.nz:])
        jac = central_difference_jacobians(dae, x, z, {})
        return jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx), dae, z

    def critical(self, A, dae, z):
        from ibr_cycles.certification.symmetry import rotation_generator
        from ibr_cycles.certification.transverse import transverse_operator

        r_x, _ = rotation_generator(dae, z)
        w = None
        if self.variant == "SP33":
            from ibr_cycles.certification.symmetry import frequency_partner
            w = frequency_partner(dae).w
        ev = np.linalg.eigvals(transverse_operator(A, r_x, w).a_perp)
        lam_t = complex(ev[np.argmax(ev.real)])
        if lam_t.imag < 0:
            lam_t = lam_t.conjugate()
        vals, vl, vr = sla.eig(A, left=True, right=True)
        j = int(np.argmin(np.abs(vals - lam_t)))
        others = ev[(np.abs(ev - lam_t) > 1e-8) & (np.abs(ev - lam_t.conjugate()) > 1e-8)]
        gap = float(lam_t.real - others.real.max()) if others.size else float("inf")
        cond = float(np.linalg.norm(vl[:, j]) * np.linalg.norm(vr[:, j]) / abs(np.vdot(vl[:, j], vr[:, j])))
        return complex(vals[j]), vl[:, j], vr[:, j], gap, float(abs(vals[j] - lam_t)), cond

    @staticmethod
    def em_top(A):
        """Rightmost EM-band eigen-triplet of A (post-hoc exploratory target)."""
        vals, vl, vr = sla.eig(A, left=True, right=True)
        f = np.abs(vals.imag) / (2 * np.pi)
        band = [j for j in range(vals.size) if R.MODE_BAND[0] <= f[j] <= R.MODE_BAND[1] and vals[j].imag > 0]
        if not band:
            return None
        j = max(band, key=lambda i: vals[i].real)
        return complex(vals[j]), vl[:, j], vr[:, j]

    def derivatives(self, branches, *, h_rel=1e-3):
        """d lambda_c / d gamma_e (frozen and total) at gamma = 1 for every branch e; also (post hoc) for the
        rightmost EM-band mode, from the same perturbed matrices."""
        w0 = self.w0()
        R0 = float(np.abs(self.residual(w0, None)).max())
        A0, dae0, z0 = self.A_of(w0, None)
        lam, vl, vr, gap, match_err, cond = self.critical(A0, dae0, z0)
        denom = np.vdot(vl, vr)
        emt = self.em_top(A0)
        if emt is not None:
            lam_em, vl_em, vr_em = emt
            den_em = np.vdot(vl_em, vr_em)
        Rw = self.jac_R(w0, None)
        lu = sla.lu_factor(Rw)
        out = []
        for e in branches:
            h = h_rel
            ap, am = {"idx": e, "value": 1.0 + h}, {"idx": e, "value": 1.0 - h}
            Ra = (self.residual(w0, ap) - self.residual(w0, am)) / (2 * h)
            dw = -sla.lu_solve(lu, Ra)
            App, _, _ = self.A_of(w0, ap)
            Amm, _, _ = self.A_of(w0, am)
            Atp, _, _ = self.A_of(w0 + h * dw, ap)
            Atm, _, _ = self.A_of(w0 - h * dw, am)
            d_fro = complex(np.vdot(vl, ((App - Amm) / (2 * h)) @ vr) / denom)
            d_tot = complex(np.vdot(vl, ((Atp - Atm) / (2 * h)) @ vr) / denom)
            rec = {"e": int(e), "d_frozen": d_fro.real, "d_frozen_im": d_fro.imag, "d_total": d_tot.real, "d_total_im": d_tot.imag,
                   "dw_norm": float(np.abs(dw).max())}
            if emt is not None:
                rec["em_d_frozen"] = float((np.vdot(vl_em, ((App - Amm) / (2 * h)) @ vr_em) / den_em).real)
                rec["em_d_total"] = float((np.vdot(vl_em, ((Atp - Atm) / (2 * h)) @ vr_em) / den_em).real)
            out.append(rec)
        return {"lam": [lam.real, lam.imag], "gap2": gap, "R0": R0, "match_err": match_err, "eig_cond": cond,
                "lam_em": [lam_em.real, lam_em.imag] if emt is not None else None, "items": out}


def finite_tracked(members, *, g, k, variant, e, gamma, lam_ref, lam_em=None, draw=None):
    """Finite re-equilibrated truth (build68 with the scaled network; independent of the engine)."""
    try:
        case = R.build68(members, variant=variant, g=g, k=k, draw=draw, network=R.network_with(scale={int(e): float(gamma)}))
    except Exception as ex:  # noqa: BLE001
        return {"ok": False, "error": repr(ex)[:200]}
    from ibr_cycles.certification.symmetry import rotation_generator
    from ibr_cycles.certification.transverse import transverse_operator

    A, _, jac = R.matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    ev = np.linalg.eigvals(transverse_operator(A, r_x, R.quotient_w(case)).a_perp)
    j = int(np.argmin(np.abs(ev - lam_ref)))
    v = case.dae.power_flow.voltages
    vm = np.abs(v)
    em = {}
    if lam_em is not None:
        vals = np.linalg.eigvals(A)
        f = np.abs(vals.imag) / (2 * np.pi)
        band = vals[(f >= R.MODE_BAND[0]) & (f <= R.MODE_BAND[1]) & (vals.imag > 0)]
        if band.size:
            jj = int(np.argmin(np.abs(band - lam_em)))
            em = {"em_tracked_re": float(band[jj].real), "em_tracked_hz": float(abs(band[jj].imag) / (2 * np.pi)), "em_top_re": float(band.real.max())}
    return {"ok": True, "alpha": float(ev.real.max()), "tracked_re": float(ev[j].real), "tracked_hz": float(abs(ev[j].imag) / (2 * np.pi)), **em,
            "vmin": float(vm.min()), "vmax": float(vm.max()), "n_v_out": int(((vm < 0.8) | (vm > 1.2)).sum()),
            "eq_residual": float(max(case.equilibrium.norm_f, case.equilibrium.norm_g))}
