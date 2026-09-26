# ruff: noqa: E501  -- sympy equation strings kept on one line
"""SG68: exact ANDES transcription of the Singh & Pal (2013) v3.3 68-bus machine with its excitation, PSS,
PST tg-model-1 governor and rotor damping (src/ibr_cycles/models/ieee68_devices.py + code/_r68.Gov68).

Three classes share one equation set and differ only in the excitation:
    SG68D  DC4B (G1-G8, G10-G12)     SG68S  ST1A (G9)     SG68M  manual Efd (G13-G16, no PSS)
Each device is on its OWN base (Singh & Pal mac_con MVA), its bus withdrawal scaled by w = Sn/100.
Machine frame as in the report: V exp(-j delta) = Vq + j Vd, injected current (Iq + j Id) exp(j delta).
The bus voltage is polar, V = v exp(j a). The internal state 'omega' (slip) is called 'slip' here.

Registered at RUNTIME only (code/andes68/run68.py) as ('cdw68_models', ['SG68D', 'SG68S', 'SG68M']) with a
private pycode folder; the shared ANDES copy of .venv/xtool-andes-gfl is not modified.
"""

from andes.core import ConstService, ExtAlgeb, ExtService, IdxParam, Model, ModelData, NumParam, State

W_B = "(2 * pi * 60)"
VQ = "(v * cos(a - delta))"
VD = "(v * sin(a - delta))"
EMFQ = f"(eq1 * kd1 + psi1d * kd2 - {VQ})"
EMFD = f"(ed1 * kq1 - psi2q * kq2 - {VD} + edc)"
DEN = "(ra**2 + xd2**2)"
IQ = f"(({EMFQ} * ra + {EMFD} * xd2) / {DEN})"
ID = f"(({EMFD} * ra - {EMFQ} * xd2) / {DEN})"
TE = f"(ed1 * {ID} * kq1 + eq1 * {IQ} * kd1 - {ID} * {IQ} * (xd2 - xq2) + psi1d * {IQ} * kd2 - psi2q * {ID} * kq2)"
PT = f"({VQ} * {IQ} + {VD} * {ID})"
QT = f"({VD} * {IQ} - {VQ} * {ID})"
# PSS: washout output and three lead-lag stages (identity stages T1 == T2 keep a decoupled state)
S0 = "(K * slip - xw)"
S1 = f"((T11 / T12) * {S0} + (1 - T11 / T12) * z1)"
S2 = f"((T21 / T22) * {S1} + (1 - T21 / T22) * z2)"
VSS = f"((T31 / T32) * {S2} + (1 - T31 / T32) * z3)"
# governor (PST tg model 1)
A1 = "(T3G / TCG)"
TM = f"(tg3 + (T4G / T5G) * (tg2 + {A1} * tg1))"


class SG68Data(ModelData):
    def __init__(self, exciter):
        super().__init__()
        self.bus = IdxParam(model="ACNode", mandatory=True, status_parent=True, info="bus")
        self.gen = IdxParam(model="StaticGen", mandatory=True, replaces=True, info="static gen")
        self.w = NumParam(default=1.0, non_zero=True, info="weight Sn/100")
        for nm, dv in (("xl", 0.0), ("ra", 0.0), ("xd", 1.0), ("xd1", 0.3), ("xd2", 0.2), ("Td10", 5.0), ("Td20", 0.05),
                       ("xq", 1.0), ("xq1", 0.3), ("xq2", 0.2), ("Tq10", 1.0), ("Tq20", 0.05), ("H", 3.0), ("D", 0.0), ("Tcd", 0.01)):
            setattr(self, nm, NumParam(default=dv, info=nm))
        for nm, dv in (("KG", 0.0), ("TSG", 0.1), ("TCG", 0.5), ("T3G", 0.0), ("T4G", 1.25), ("T5G", 5.0)):
            setattr(self, nm, NumParam(default=dv, info=f"governor {nm}"))
        if exciter == "DC4B":
            for nm, dv in (("TR", 0.01), ("KA", 1.0), ("GS", 1.0), ("TA", 0.02), ("KE", 1.0), ("TE", 0.785), ("AEX", 0.0), ("BEX", 0.0),
                           ("KF", 0.03), ("TF", 1.0), ("KP", 200.0), ("KI", 50.0), ("KD", 50.0), ("TD", 0.01)):
                setattr(self, nm, NumParam(default=dv, info=f"DC4B {nm}"))
        elif exciter == "ST1A":
            for nm, dv in (("TR", 0.01), ("KA", 200.0), ("GS", 1.0)):
                setattr(self, nm, NumParam(default=dv, info=f"ST1A {nm}"))
        if exciter in ("DC4B", "ST1A"):
            for nm, dv in (("K", 20.0), ("TW", 15.0), ("T11", 0.15), ("T12", 0.04), ("T21", 0.15), ("T22", 0.04), ("T31", 0.15), ("T32", 0.04)):
                setattr(self, nm, NumParam(default=dv, info=f"PSS {nm}"))


class SG68Model(Model):
    def __init__(self, system, config, exciter):
        Model.__init__(self, system, config)
        self.flags.update({"tds": True, "nr_iter": False})
        self.a = ExtAlgeb(model="Bus", src="a", indexer=self.bus, e_str=f"-ue * w * {PT}")
        self.v = ExtAlgeb(model="Bus", src="v", indexer=self.bus, e_str=f"-ue * w * {QT}")
        self.kd1 = ConstService(v_str="(xd2 - xl) / (xd1 - xl)")
        self.kd2 = ConstService(v_str="(xd1 - xd2) / (xd1 - xl)")
        self.kq1 = ConstService(v_str="(xq2 - xl) / (xq1 - xl)")
        self.kq2 = ConstService(v_str="(xq1 - xq2) / (xq1 - xl)")
        # initialization: Ieee68Machine.initialize, exactly
        self.p0 = ExtService(model="StaticGen", src="p", indexer=self.gen)
        self.q0 = ExtService(model="StaticGen", src="q", indexer=self.gen)
        self._V = ConstService(v_str="v * exp(1j * a)", vtype=complex)
        self._I = ConstService(v_str="conj((p0 + 1j * q0) / w / _V)", vtype=complex)
        self._E = ConstService(v_str="_V + (ra + 1j * xq) * _I", vtype=complex)
        self._dc = ConstService(v_str="log(_E / abs(_E))", vtype=complex)
        self.delta0 = ConstService(v_str="im(_dc)")
        self._vr = ConstService(v_str="_V * exp(-_dc)", vtype=complex)
        self._ir = ConstService(v_str="_I * exp(-_dc)", vtype=complex)
        self.vq0 = ConstService(v_str="re(_vr)")
        self.vd0 = ConstService(v_str="im(_vr)")
        self.iq0 = ConstService(v_str="re(_ir)")
        self.id0 = ConstService(v_str="im(_ir)")
        self.ed10 = ConstService(v_str="-(xq - xq1) * iq0")
        self.psi2q0 = ConstService(v_str="-ed10 + (xq1 - xl) * iq0")
        self.edc0 = ConstService(v_str="iq0 * (xd2 - xq2)")
        self.eq10 = ConstService(v_str="vq0 - xd1 * id0 + ra * iq0")
        self.psi1d0 = ConstService(v_str="eq10 + (xd1 - xl) * id0")
        self.efd0 = ConstService(v_str="eq10 - (xd - xd1) * id0")
        self.pm0 = ConstService(v_str="vq0 * iq0 + vd0 * id0 + ra * (iq0**2 + id0**2)")
        self.delta = State(v_str="delta0", e_str=f"ue * {W_B} * slip")
        self.slip = State(v_str="0", e_str=f"ue * ({TM} - {TE} - D * slip) / (2 * H)")
        if exciter == "DC4B":
            efd = "efd"
        elif exciter == "ST1A":
            efd = f"(GS * KA * (vref + {VSS} - vr))"
        else:
            efd = "efd0"
        self.eq1 = State(v_str="eq10", e_str=f"ue * ({efd} - eq1 + (xd - xd1) * ({ID} + (xd1 - xd2) / (xd1 - xl)**2 * (psi1d - (xd1 - xl) * {ID} - eq1))) / Td10")
        self.ed1 = State(v_str="ed10", e_str=f"ue * (-ed1 + (xq - xq1) * (-{IQ} + (xq1 - xq2) / (xq1 - xl)**2 * ((xq1 - xl) * {IQ} - ed1 - psi2q))) / Tq10")
        self.psi1d = State(v_str="psi1d0", e_str=f"ue * (eq1 + (xd1 - xl) * {ID} - psi1d) / Td20")
        self.psi2q = State(v_str="psi2q0", e_str=f"ue * (-ed1 + (xq1 - xl) * {IQ} - psi2q) / Tq20")
        self.edc = State(v_str="edc0", e_str=f"ue * ({IQ} * (xd2 - xq2) - edc) / Tcd")
        if exciter == "DC4B":
            self.va0 = ConstService(v_str="KE * efd0 + efd0 * AEX * exp(BEX * efd0)")
            self.vref = ConstService(v_str="v")
            u = f"(vref + {VSS} - vr - (KF / TF) * (efd - vf))"
            vpid = f"(KP * {u} + xi + (KD / TD) * ({u} - xdd))"
            self.vr = State(v_str="v", e_str="ue * (v - vr) / TR")
            self.vf = State(v_str="efd0", e_str="ue * (efd - vf) / TF")
            self.va = State(v_str="va0", e_str=f"ue * (GS * KA * {vpid} - va) / TA")
            self.efd = State(v_str="efd0", e_str="ue * (va - (KE * efd + efd * AEX * exp(BEX * efd))) / TE")
            self.xi = State(v_str="va0 / (GS * KA)", e_str=f"ue * KI * {u}")
            self.xdd = State(v_str="0", e_str=f"ue * ({u} - xdd) / TD")
        elif exciter == "ST1A":
            self.vref = ConstService(v_str="v + efd0 / (GS * KA)")
            self.vr = State(v_str="v", e_str="ue * (v - vr) / TR")
        if exciter in ("DC4B", "ST1A"):
            self.xw = State(v_str="0", e_str=f"ue * (K * slip - xw) / TW")
            self.z1 = State(v_str="0", e_str=f"ue * ({S0} - z1) / T12")
            self.z2 = State(v_str="0", e_str=f"ue * ({S1} - z2) / T22")
            self.z3 = State(v_str="0", e_str=f"ue * ({S2} - z3) / T32")
        self.tg1 = State(v_str="pm0", e_str="ue * (pm0 - KG * slip - tg1) / TSG")
        self.tg2 = State(v_str=f"(1 - {A1}) * pm0", e_str=f"ue * ((1 - {A1}) * tg1 - tg2) / TCG")
        self.tg3 = State(v_str="(1 - T4G / T5G) * pm0", e_str=f"ue * ((tg2 + {A1} * tg1) * (1 - T4G / T5G) - tg3) / T5G")


class SG68D(SG68Data, SG68Model):
    """Singh & Pal machine + DC4B + speed PSS + PST governor."""

    def __init__(self, system, config):
        SG68Data.__init__(self, "DC4B")
        SG68Model.__init__(self, system, config, "DC4B")


class SG68S(SG68Data, SG68Model):
    """Singh & Pal machine + ST1A + speed PSS + PST governor."""

    def __init__(self, system, config):
        SG68Data.__init__(self, "ST1A")
        SG68Model.__init__(self, system, config, "ST1A")


class SG68M(SG68Data, SG68Model):
    """Singh & Pal machine + manual Efd + PST governor (area equivalents)."""

    def __init__(self, system, config):
        SG68Data.__init__(self, "MANUAL")
        SG68Model.__init__(self, system, config, "MANUAL")
