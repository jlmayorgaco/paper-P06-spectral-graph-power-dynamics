# ruff: noqa: E501  -- sympy equation strings and parameter tables kept on one line
"""Custom ANDES models for the PCV Phase 8 reproduction (spec commit 3e847a4f).

Exact transcriptions of the internal device equations
(src/ibr_cycles/models/ieee39_devices.py): each device is on its OWN base and its
bus withdrawal is scaled by the rating weight ``w = Sn / 100``. No ANDES per-unit
conversion flag is used. The bus voltage is polar, V = v exp(j a).

This file is copied into the isolated ANDES copy as ``andes/models/pcv_models.py``
and registered as ``('pcv_models', ['SG2AX', 'GFL11'])``. It is never imported by
the internal code base.
"""

from andes.core import (
    ConstService,
    ExtAlgeb,
    ExtService,
    IdxParam,
    Model,
    ModelData,
    NumParam,
    State,
)

W_B = "(2 * pi * 60)"  # OMEGA_B of the internal model

# ----------------------------------------------------------------- SG2AX ----
SG_VD = "(v * sin(delta - a))"
SG_VQ = "(v * cos(delta - a))"
SG_DET = "(ra**2 + xd1 * xq1)"
SG_ID = f"((ra * (ed1 - {SG_VD}) + xq1 * (eq1 - {SG_VQ})) / {SG_DET})"
SG_IQ = f"((-xd1 * (ed1 - {SG_VD}) + ra * (eq1 - {SG_VQ})) / {SG_DET})"
SG_PT = f"({SG_VD} * {SG_ID} + {SG_VQ} * {SG_IQ})"
SG_QT = f"({SG_VQ} * {SG_ID} - {SG_VD} * {SG_IQ})"
SG_PE = f"({SG_PT} + ra * ({SG_ID}**2 + {SG_IQ}**2))"


class SG2AXData(ModelData):
    def __init__(self):
        super().__init__()
        self.bus = IdxParam(
            model="ACNode", mandatory=True, status_parent=True, info="bus"
        )
        self.gen = IdxParam(
            model="StaticGen", mandatory=True, replaces=True, info="static gen"
        )
        self.w = NumParam(default=1.0, non_zero=True, info="rating weight Sn/100")
        self.ra = NumParam(default=0.0, info="armature resistance (machine base)")
        self.xd = NumParam(default=1.0, info="xd")
        self.xq = NumParam(default=1.0, info="xq")
        self.xd1 = NumParam(default=0.3, non_zero=True, info="xd'")
        self.xq1 = NumParam(default=0.3, non_zero=True, info="xq'")
        self.Td10 = NumParam(default=5.0, non_zero=True, info="Td0'")
        self.Tq10 = NumParam(default=1.0, non_zero=True, info="Tq0'")
        self.M = NumParam(default=6.0, non_zero=True, info="2H (machine base)")
        self.D = NumParam(default=0.0, info="damping (machine base)")
        self.KA = NumParam(default=1.0, non_zero=True, info="AVR gain (policy-scaled)")
        self.TE = NumParam(
            default=1.0, non_zero=True, info="AVR time constant (policy-scaled)"
        )
        self.KS = NumParam(default=0.0, info="PSS gain")
        self.T4 = NumParam(default=1.0, non_zero=True, info="PSS lag")
        self.T5 = NumParam(default=1.0, info="PSS washout gain time constant")
        self.T6 = NumParam(default=1.0, non_zero=True, info="PSS washout lag")


class SG2AXModel(Model):
    def __init__(self, system, config):
        Model.__init__(self, system, config)
        self.flags.update({"tds": True, "nr_iter": False})

        self.a = ExtAlgeb(
            model="Bus", src="a", indexer=self.bus, e_str=f"-ue * w * {SG_PT}"
        )
        self.v = ExtAlgeb(
            model="Bus", src="v", indexer=self.bus, e_str=f"-ue * w * {SG_QT}"
        )

        # initialization (spec 4.1)
        self.p0 = ExtService(model="StaticGen", src="p", indexer=self.gen)
        self.q0 = ExtService(model="StaticGen", src="q", indexer=self.gen)
        self._V = ConstService(v_str="v * exp(1j * a)", vtype=complex)
        self._I = ConstService(v_str="conj((p0 + 1j * q0) / w / _V)", vtype=complex)
        self._E = ConstService(v_str="_V + (ra + 1j * xq) * _I", vtype=complex)
        self._dc = ConstService(v_str="log(_E / abs(_E))", vtype=complex)
        self.delta0 = ConstService(v_str="im(_dc)")
        self._vdq = ConstService(v_str="_V * exp(1j * 0.5 * pi - _dc)", vtype=complex)
        self._idq = ConstService(v_str="_I * exp(1j * 0.5 * pi - _dc)", vtype=complex)
        self.vd0 = ConstService(v_str="re(_vdq)")
        self.vq0 = ConstService(v_str="im(_vdq)")
        self.id0 = ConstService(v_str="re(_idq)")
        self.iq0 = ConstService(v_str="im(_idq)")
        self.eq10 = ConstService(v_str="vq0 + ra * iq0 + xd1 * id0")
        self.ed10 = ConstService(v_str="vd0 + ra * id0 - xq1 * iq0")
        self.efd0 = ConstService(v_str="eq10 + (xd - xd1) * id0")
        self.pm = ConstService(v_str="vd0 * id0 + vq0 * iq0 + ra * (id0**2 + iq0**2)")
        self.vref = ConstService(v_str="v + efd0 / KA")

        self.delta = State(v_str="delta0", e_str=f"ue * {W_B} * (omega - 1)")
        self.omega = State(
            v_str="u", e_str=f"ue * (pm - {SG_PE} - D * (omega - 1)) / M"
        )
        self.eq1 = State(
            v_str="eq10", e_str=f"ue * (efd - eq1 - (xd - xd1) * {SG_ID}) / Td10"
        )
        self.ed1 = State(
            v_str="ed10", e_str=f"ue * (-ed1 + (xq - xq1) * {SG_IQ}) / Tq10"
        )
        self.efd = State(
            v_str="efd0", e_str="ue * (KA * (vref + KS * pl - v) - efd) / TE"
        )
        self.pw = State(v_str="pm", e_str=f"ue * ({SG_PE} - pw) / T6")
        self.pl = State(v_str="0", e_str=f"ue * (T5 * ({SG_PE} - pw) / T6 - pl) / T4")


class SG2AX(SG2AXData, SG2AXModel):
    """Two-axis machine + first-order AVR + washout-lag power PSS (7 states)."""

    def __init__(self, system, config):
        SG2AXData.__init__(self)
        SG2AXModel.__init__(self, system, config)


# ----------------------------------------------------------------- GFL11 ----
GF_VD = "(v * cos(a - theta))"
GF_VQ = "(v * sin(a - theta))"
GF_P = f"({GF_VD} * i_d + {GF_VQ} * i_q)"
GF_Q = f"({GF_VQ} * i_d - {GF_VD} * i_q)"
GF_ERR = "(g * (v_ref - v))"
GF_QCMD = f"(kp_v * {GF_ERR} + x_v)"
GF_IDREF = "(kp_p * (p_ref - p_f) + x_p)"
GF_IQREF = f"(-(kp_q * ({GF_QCMD} - q_f) + x_q))"
GF_ED = f"({GF_VD} + kp_i * ({GF_IDREF} - i_d) + x_id - xf * i_q)"
GF_EQ = f"({GF_VQ} + kp_i * ({GF_IQREF} - i_q) + x_iq + xf * i_d)"


class GFL11Data(ModelData):
    def __init__(self):
        super().__init__()
        self.bus = IdxParam(
            model="ACNode", mandatory=True, status_parent=True, info="bus"
        )
        self.gen = IdxParam(
            model="StaticGen", mandatory=True, replaces=True, info="static gen"
        )
        self.w = NumParam(default=1.0, non_zero=True, info="rating weight Sn/100")
        self.kp_pll = NumParam(default=53.0, info="PLL kp")
        self.ki_pll = NumParam(default=1400.0, info="PLL ki")
        self.tau_p = NumParam(default=0.03, non_zero=True, info="power filter")
        self.kp_p = NumParam(default=0.20, info="P loop kp")
        self.ki_p = NumParam(default=8.0, info="P loop ki")
        self.kp_q = NumParam(default=0.20, info="Q loop kp")
        self.ki_q = NumParam(default=8.0, info="Q loop ki")
        self.kp_i = NumParam(default=0.25, info="current loop kp")
        self.ki_i = NumParam(default=6.0, info="current loop ki")
        self.xf = NumParam(default=0.15, non_zero=True, info="filter reactance")
        self.rf = NumParam(default=0.01, info="filter resistance")
        self.kp_v = NumParam(default=2.0, info="Q/V kp")
        self.ki_v = NumParam(default=20.0, info="Q/V ki")
        self.g = NumParam(default=1.0, info="policy Q/V gain")
        self.leak = NumParam(default=0.05, info="Q/V integrator leak (rad/s)")


class GFL11Model(Model):
    def __init__(self, system, config):
        Model.__init__(self, system, config)
        self.flags.update({"tds": True, "nr_iter": False})

        self.a = ExtAlgeb(
            model="Bus", src="a", indexer=self.bus, e_str=f"-ue * w * {GF_P}"
        )
        self.v = ExtAlgeb(
            model="Bus", src="v", indexer=self.bus, e_str=f"-ue * w * {GF_Q}"
        )

        # initialization (spec 4.2)
        self.p0 = ExtService(model="StaticGen", src="p", indexer=self.gen)
        self.q0 = ExtService(model="StaticGen", src="q", indexer=self.gen)
        self._V = ConstService(v_str="v * exp(1j * a)", vtype=complex)
        self._I = ConstService(v_str="conj((p0 + 1j * q0) / w / _V)", vtype=complex)
        self._idq = ConstService(v_str="_I * exp(-1j * a)", vtype=complex)
        self.id0 = ConstService(v_str="re(_idq)")
        self.iq0 = ConstService(v_str="im(_idq)")
        self.p_ref = ConstService(v_str="v * id0")
        self.q_ref = ConstService(v_str="-v * iq0")
        self.v_ref = ConstService(v_str="v")

        self.theta = State(v_str="a", e_str=f"ue * (kp_pll * {GF_VQ} + x_pll)")
        self.x_pll = State(v_str="0", e_str=f"ue * ki_pll * {GF_VQ}")
        self.p_f = State(v_str="p_ref", e_str=f"ue * ({GF_P} - p_f) / tau_p")
        self.q_f = State(v_str="q_ref", e_str=f"ue * ({GF_Q} - q_f) / tau_p")
        self.x_p = State(v_str="id0", e_str="ue * ki_p * (p_ref - p_f)")
        self.x_q = State(v_str="-iq0", e_str=f"ue * ki_q * ({GF_QCMD} - q_f)")
        self.i_d = State(
            v_str="id0",
            e_str=f"ue * ({W_B} / xf) * ({GF_ED} - {GF_VD} - rf * i_d + xf * i_q)",
        )
        self.i_q = State(
            v_str="iq0",
            e_str=f"ue * ({W_B} / xf) * ({GF_EQ} - {GF_VQ} - rf * i_q - xf * i_d)",
        )
        self.x_id = State(v_str="rf * id0", e_str=f"ue * ki_i * ({GF_IDREF} - i_d)")
        self.x_iq = State(v_str="rf * iq0", e_str=f"ue * ki_i * ({GF_IQREF} - i_q)")
        self.x_v = State(
            v_str="q_ref", e_str=f"ue * (ki_v * {GF_ERR} - leak * (x_v - q_ref))"
        )


class GFL11(GFL11Data, GFL11Model):
    """Grid-following converter with Q/V regulator and leak (11 states)."""

    def __init__(self, system, config):
        GFL11Data.__init__(self)
        GFL11Model.__init__(self, system, config)
