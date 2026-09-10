"""Documented machine, exciter and PSS models of the IEEE 68-bus / NETS-NYPS benchmark.

Source: Singh & Pal (2013), "Report on the 68-bus, 16-machine, 5-area system",
IEEE PES Task Force on Benchmark Systems for Stability Controls, v3.3, Section
II.A, equations (1)-(10), and the initialization in its Init_MultiMachine.m.
They are implemented as printed:

    machine   sub-transient model with four rotor coils (E'q, E'd, psi1d, psi2q)
              plus the dummy coil E'dc (Tc = 0.01 s) for sub-transient saliency;
              equations (1)-(6) and the stator equation
    DC4B      equation (7): input filter, PID regulator with rate feedback,
              regulator lag, DC exciter with exponential saturation
    ST1A      equation (8): input filter and static gain
    MANUAL    field voltage held at its initial value (G13-G16)
    PSS       equation (9): slip input, washout and three lead-lags

The machine frame is the report's, with the q axis real and the d axis
imaginary: V exp(-j delta) = Vq + j Vd. The stator current (Iq + j Id)
exp(j delta) is the current injected into the network. Every quantity is per
unit on the machine's MVA base (the report's BM = MVA/100 conversion is the
device ``weight``), so the device exposes the same interface as the IEEE-39
devices.

Controller limits are not represented, since this is a small-signal model.
``initialize`` records how far the equilibrium sits from each documented limit
(``limit_margin``). The slip S = omega - 1 is carried under the label
``omega_sg<bus>``, so that code looking for machine-speed states finds it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

OMEGA_B = 2.0 * math.pi * 60.0
MACHINE_STATES = ("delta", "omega", "eq1", "ed1", "psi1d", "psi2q", "edc")
DC4B_STATES = ("vr", "vf", "va", "efd", "pid_i", "pid_d")


@dataclass(frozen=True)
class Ieee68MachineParameters:
    """One documented generator with its excitation and stabilizer (machine base)."""

    xl: float
    ra: float
    xd: float
    xd1: float
    xd2: float
    td10: float
    td20: float
    xq: float
    xq1: float
    xq2: float
    tq10: float
    tq20: float
    h: float
    d: float
    exciter: str  # "DC4B" | "ST1A" | "MANUAL"
    exc: tuple[tuple[str, float], ...] = ()
    pss: tuple[tuple[str, float], ...] = ()
    tc: float = 0.01
    #: excitation coordinate of the replacement study: scales the regulator gain
    #: (DC4B Ka, ST1A KA). 1 is the documented case.
    gain_scale: float = 1.0
    # setpoints, filled at initialization
    tm: float = 0.0
    vref: float = 1.0
    efd0: float = 0.0

    @classmethod
    def from_rows(cls, machine: dict, avr: dict, pss: dict) -> Ieee68MachineParameters:
        exc_type = avr.get("type", "MANUAL")
        exc: dict[str, float] = {}
        if exc_type == "DC4B":
            e1, se1, e2, se2 = avr["E1"], avr["SE1"], avr["E2"], avr["SE2"]
            bex = math.log(se1 / se2) / (e1 - e2)
            aex = se1 * math.exp(-bex * e1)
            # Init_MultiMachine.m: Kd <- column 17, Ki <- column 18 (both 50 in v3.3)
            exc = {
                "tr": avr["TR"],
                "ka": avr["KA"],
                "ta": avr["TA"],
                "ke": avr["KE"],
                "te": avr["TE"],
                "aex": aex,
                "bex": bex,
                "kf": avr["KF"],
                "tf": avr["TF"],
                "kp": avr["KP"],
                "ki": avr["KD"],
                "kd": avr["KI"],
                "td": avr["TD"],
                "vrmax": avr["VRmax"],
                "vrmin": avr["VRmin"],
            }
        elif exc_type == "ST1A":
            exc = {
                "tr": avr["TR"],
                "ka": avr["KA"],
                "vrmax": avr["VRmax"],
                "vrmin": avr["VRmin"],
            }
        stab: dict[str, float] = {}
        if pss.get("type") == "SPEED3":
            stab = {
                k.lower(): float(pss[k])
                for k in (
                    "K",
                    "TW",
                    "T11",
                    "T12",
                    "T21",
                    "T22",
                    "T31",
                    "T32",
                    "VSmax",
                    "VSmin",
                )
            }
        return cls(
            xl=machine["xl"],
            ra=machine["ra"],
            xd=machine["xd"],
            xd1=machine["xd1"],
            xd2=machine["xd2"],
            td10=machine["Td10"],
            td20=machine["Td20"],
            xq=machine["xq"],
            xq1=machine["xq1"],
            xq2=machine["xq2"],
            tq10=machine["Tq10"],
            tq20=machine["Tq20"],
            h=machine["H"],
            d=machine["d0"],
            exciter=exc_type,
            exc=tuple(sorted((k, float(v)) for k, v in exc.items())),
            pss=tuple(sorted(stab.items())),
        )

    @property
    def e(self) -> dict[str, float]:
        return dict(self.exc)

    @property
    def p(self) -> dict[str, float]:
        return dict(self.pss)

    @property
    def lead_lags(self) -> tuple[tuple[float, float], ...]:
        """Stabilizer stages actually carrying dynamics (T1 == T2 is an identity)."""

        p = self.p
        if not p:
            return ()
        stages = ((p["t11"], p["t12"]), (p["t21"], p["t22"]), (p["t31"], p["t32"]))
        return tuple(s for s in stages if s[0] != s[1])


@dataclass
class Ieee68Machine:
    """Documented 68-bus generator as a bus-connected device."""

    bus: int
    parameters: Ieee68MachineParameters
    weight: float = 1.0
    limit_margin: float = float("inf")

    # ----------------------------------------------------------- structure --
    @property
    def _exc_states(self) -> tuple[str, ...]:
        if self.parameters.exciter == "DC4B":
            return DC4B_STATES
        if self.parameters.exciter == "ST1A":
            return ("vr",)
        return ()

    @property
    def _pss_states(self) -> tuple[str, ...]:
        if not self.parameters.pss:
            return ()
        return ("pss_w",) + tuple(
            f"pss_l{i + 1}" for i in range(len(self.parameters.lead_lags))
        )

    @property
    def labels(self) -> tuple[str, ...]:
        names = MACHINE_STATES + self._exc_states + self._pss_states
        return tuple(f"{name}_sg{self.bus}" for name in names)

    @property
    def n_states(self) -> int:
        return len(MACHINE_STATES) + len(self._exc_states) + len(self._pss_states)

    # ------------------------------------------------------------- machine --
    def _stator(self, x, v: complex):
        p = self.parameters
        delta, eq1, ed1, psi1d, psi2q, edc = x[0], x[2], x[3], x[4], x[5], x[6]
        rotated = v * complex(math.cos(delta), -math.sin(delta))
        vq, vd = rotated.real, rotated.imag
        kd1 = (p.xd2 - p.xl) / (p.xd1 - p.xl)
        kd2 = (p.xd1 - p.xd2) / (p.xd1 - p.xl)
        kq1 = (p.xq2 - p.xl) / (p.xq1 - p.xl)
        kq2 = (p.xq1 - p.xq2) / (p.xq1 - p.xl)
        emf = complex(eq1 * kd1 + psi1d * kd2 - vq, ed1 * kq1 - psi2q * kq2 - vd + edc)
        current = emf / complex(p.ra, p.xd2)
        iq, id_ = current.real, current.imag
        te = (
            ed1 * id_ * kq1
            + eq1 * iq * kd1
            - id_ * iq * (p.xd2 - p.xq2)
            + psi1d * iq * kd2
            - psi2q * id_ * kq2
        )
        return iq, id_, te

    def _efd(self, x) -> float:
        p = self.parameters
        n = len(MACHINE_STATES)
        if p.exciter == "DC4B":
            return float(x[n + 3])
        if p.exciter == "ST1A":
            e = p.e
            return p.gain_scale * e["ka"] * (p.vref + self._vss(x) - float(x[n]))
        return p.efd0

    def _vss(self, x) -> float:
        p = self.parameters
        if not p.pss:
            return 0.0
        s = p.p
        start = len(MACHINE_STATES) + len(self._exc_states)
        slip = float(x[1])
        signal = s["k"] * slip - float(x[start])  # washout output
        for i, (t1, t2) in enumerate(p.lead_lags):
            z = float(x[start + 1 + i])
            signal = (t1 / t2) * signal + (1.0 - t1 / t2) * z
        return signal

    def derivatives(self, x: NDArray[np.float64], v: complex) -> NDArray[np.float64]:
        p = self.parameters
        iq, id_, te = self._stator(x, v)
        slip, eq1, ed1, psi1d, psi2q, edc = (float(x[i]) for i in range(1, 7))
        efd = self._efd(x)
        out = np.empty(self.n_states)
        out[0] = OMEGA_B * slip
        out[1] = (p.tm - te - p.d * slip) / (2.0 * p.h)
        out[2] = (
            efd
            - eq1
            + (p.xd - p.xd1)
            * (
                id_
                + (p.xd1 - p.xd2)
                / (p.xd1 - p.xl) ** 2
                * (psi1d - (p.xd1 - p.xl) * id_ - eq1)
            )
        ) / p.td10
        out[3] = (
            -ed1
            + (p.xq - p.xq1)
            * (
                -iq
                + (p.xq1 - p.xq2)
                / (p.xq1 - p.xl) ** 2
                * ((p.xq1 - p.xl) * iq - ed1 - psi2q)
            )
        ) / p.tq10
        out[4] = (eq1 + (p.xd1 - p.xl) * id_ - psi1d) / p.td20
        out[5] = (-ed1 + (p.xq1 - p.xl) * iq - psi2q) / p.tq20
        out[6] = (iq * (p.xd2 - p.xq2) - edc) / p.tc
        n = len(MACHINE_STATES)
        terminal = abs(v)
        vss = self._vss(x)
        if p.exciter == "DC4B":
            e = p.e
            vr, vf, va, efd_s, xi, xd = (float(x[n + i]) for i in range(6))
            u = p.vref + vss - vr - (e["kf"] / e["tf"]) * (efd_s - vf)
            vpid = e["kp"] * u + xi + (e["kd"] / e["td"]) * (u - xd)
            out[n + 0] = (terminal - vr) / e["tr"]
            out[n + 1] = (efd_s - vf) / e["tf"]
            out[n + 2] = (p.gain_scale * e["ka"] * vpid - va) / e["ta"]
            out[n + 3] = (
                va - (e["ke"] * efd_s + efd_s * e["aex"] * math.exp(e["bex"] * efd_s))
            ) / e["te"]
            out[n + 4] = e["ki"] * u
            out[n + 5] = (u - xd) / e["td"]
        elif p.exciter == "ST1A":
            out[n] = (terminal - float(x[n])) / p.e["tr"]
        if p.pss:
            s = p.p
            start = n + len(self._exc_states)
            out[start] = (s["k"] * slip - float(x[start])) / s["tw"]
            signal = s["k"] * slip - float(x[start])
            for i, (t1, t2) in enumerate(p.lead_lags):
                z = float(x[start + 1 + i])
                out[start + 1 + i] = (signal - z) / t2
                signal = (t1 / t2) * signal + (1.0 - t1 / t2) * z
        return out

    def injection(self, x: NDArray[np.float64], v: complex) -> complex:
        iq, id_, _ = self._stator(x, v)
        delta = float(x[0])
        return (
            self.weight * complex(iq, id_) * complex(math.cos(delta), math.sin(delta))
        )

    def scaled(self, weight: float) -> Ieee68Machine:
        return Ieee68Machine(bus=self.bus, parameters=self.parameters, weight=weight)

    # ------------------------------------------------------- initialization --
    def initialize(
        self, v: complex, s: complex
    ) -> tuple[Ieee68Machine, NDArray[np.float64]]:
        p = self.parameters
        current = np.conj(s / self.weight / v)
        internal = v + complex(p.ra, p.xq) * current
        delta = float(np.angle(internal))
        back = complex(math.cos(delta), -math.sin(delta))
        vr_, ir_ = v * back, current * back
        vq, iq, id_ = vr_.real, ir_.real, ir_.imag
        ed1 = -(p.xq - p.xq1) * iq
        psi2q = -ed1 + (p.xq1 - p.xl) * iq
        edc = iq * (p.xd2 - p.xq2)
        eq1 = vq - p.xd1 * id_ + p.ra * iq
        psi1d = eq1 + (p.xd1 - p.xl) * id_
        efd = eq1 - (p.xd - p.xd1) * id_
        terminal = float(abs(v))
        machine = [delta, 0.0, eq1, ed1, psi1d, psi2q, edc]
        margins = []
        exc_state: list[float] = []
        vref = terminal
        if p.exciter == "DC4B":
            e = p.e
            va = e["ke"] * efd + efd * e["aex"] * math.exp(e["bex"] * efd)
            vpid = va / (p.gain_scale * e["ka"])
            exc_state = [terminal, efd, va, efd, vpid, 0.0]
            margins += [
                e["vrmax"] - va,
                va - e["vrmin"],
                e["vrmax"] / (p.gain_scale * e["ka"]) - vpid,
                vpid - e["vrmin"] / (p.gain_scale * e["ka"]),
            ]
        elif p.exciter == "ST1A":
            e = p.e
            vref = terminal + efd / (p.gain_scale * e["ka"])
            exc_state = [terminal]
            margins += [e["vrmax"] - efd, efd - e["vrmin"]]
        pss_state = [0.0] * len(self._pss_states)
        tuned = replace(p, vref=vref, efd0=efd)
        device = Ieee68Machine(
            bus=self.bus,
            parameters=tuned,
            weight=self.weight,
            limit_margin=min(margins) if margins else float("inf"),
        )
        x0 = np.array(machine + exc_state + pss_state)
        _, _, te = device._stator(x0, v)
        device.parameters = replace(tuned, tm=te)
        return device, x0
