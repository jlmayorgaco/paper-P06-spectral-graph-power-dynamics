# ruff: noqa: E501
"""CDW68 model library: IEEE 68-bus replication of the CDW contextual-weakness tests.

Frozen inputs:
- network and devices: Singh & Pal (2013) v3.3 benchmark as transcribed and validated in Gate 3
  (configs/ieee68/ieee68_network.json, src/ibr_cycles/models/ieee68_devices.py), used READ-ONLY
  through ibr_cycles.models.ieee39_case.build_dae;
- primary-frequency dynamics: inputs/pst_primary_frequency_v1.json (PST data16m governor block and
  machine ratings/damping; PST tg model 1 equations);
- converter MODEL A: the frozen TX4 grid-following converter (ConverterParameters defaults,
  leaky Q/V regulator, leak 0.05 rad/s), rated |S_gen|/0.8 (G3 rule).

Variants (docs/CDW68_PREREG_V1.md): REAL (confirmatory), NOGOV (ablation B), SP33 (ablation C).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
P68 = HERE.parent  # contextual_dynamic_weakness/ieee68_replication
CDW = P68.parent
RESEARCH = CDW.parent
for p in (RESEARCH / "src", RESEARCH / "experiments"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.dynamics.equilibrium import solve_equilibrium  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.ieee39_case import InfeasibleReplacement, ReplacementPlan, build_dae  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402
from ibr_cycles.models.ieee68_devices import Ieee68Machine  # noqa: E402

NETWORK = RESEARCH / "configs" / "ieee68" / "ieee68_network.json"
INPUTS = P68 / "inputs"
RESULTS = P68 / "results"
RAW = P68 / "raw"
LOGS = P68 / "logs"
FIGS = P68 / "figures"
DOCS = P68 / "docs"
for _d in (RESULTS, RAW, LOGS, FIGS):
    _d.mkdir(parents=True, exist_ok=True)

PRIMARY = INPUTS / "pst_primary_frequency_v1.json"
VARIANTS = ("REAL", "NOGOV", "SP33")
LEAK = 0.05  # TX4 _f7_common.LEAK (frozen)
TAU_MAT = 0.01
TAU_RES = 1e-3
MAC_MIN = 0.80
DF_MAX = 0.15
MODE_BAND = (0.1, 2.0)
MODE_RE_MIN = -1.0
STRUCT_ZERO = 1e-6
MARGIN = 2e-3  # TX4 verdict band
V68 = (3, 4, 6, 9, 11, 12)  # docs/CDW68_PREREG_V1.md section 4 (sorted by bus)
PHYSICAL = tuple(range(1, 13))


def env_threads():
    for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


# ------------------------------------------------------------------- inputs --
@lru_cache(maxsize=1)
def network_payload() -> dict:
    return json.loads(NETWORK.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def primary() -> dict:
    return json.loads(PRIMARY.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def base_network():
    return load_network(NETWORK)


def sn_sp(bus: int) -> float:
    return float(next(m["Sn"] for m in network_payload()["machines"] if int(m["bus"]) == bus))


def rating(bus: int) -> float:
    return float(next(m["rating_mva"] for m in primary()["machines"] if int(m["sp_bus"]) == bus))


def pst_damping(bus: int) -> float:
    """d_o (pu on PST rating) converted to the Singh & Pal device base (same MW per pu speed)."""
    d0 = float(next(m["pst_d0"] for m in primary()["machines"] if int(m["sp_bus"]) == bus))
    return d0 * rating(bus) / sn_sp(bus)


def label(members) -> str:
    return "+".join(map(str, sorted(members))) or "BASE"


def subsets(cands=V68):
    return [tuple(sorted(s)) for r in range(len(cands) + 1) for s in combinations(cands, r)]


def mask_of(members, cands=V68) -> int:
    return sum(1 << cands.index(b) for b in members)


def members_of(mask: int, cands=V68) -> tuple:
    return tuple(b for k, b in enumerate(cands) if mask >> k & 1)


# ------------------------------------------------------------- the governor --
@dataclass(frozen=True)
class PstTg:
    kg: float  # droop gain on the device base: pu power per pu slip (1/R * rating / Sn_SP)
    ts: float
    tc: float
    t3: float
    t4: float
    t5: float


def pst_tg(bus: int) -> PstTg:
    g = primary()["governor"]
    return PstTg(kg=g["inv_R_pu_machine_base"] * rating(bus) / sn_sp(bus), ts=g["Ts"], tc=g["Tc"], t3=g["T3"], t4=g["T4"], t5=g["T5"])


@dataclass
class Gov68:
    """Singh & Pal machine with the PST tg model-1 turbine governor (tg1, tg2, tg3)."""

    base: Ieee68Machine
    tg: PstTg
    pref: float

    @property
    def bus(self) -> int:
        return self.base.bus

    @property
    def weight(self) -> float:
        return self.base.weight

    @property
    def parameters(self):
        return self.base.parameters

    @property
    def labels(self) -> tuple[str, ...]:
        b = self.base.bus
        return (*self.base.labels, f"tg1_sg{b}", f"tg2_sg{b}", f"tg3_sg{b}")

    @property
    def n_states(self) -> int:
        return self.base.n_states + 3

    def pm(self, x) -> float:
        n0 = self.base.n_states
        tg1, tg2, tg3 = float(x[n0]), float(x[n0 + 1]), float(x[n0 + 2])
        a1 = self.tg.t3 / self.tg.tc
        return tg3 + (self.tg.t4 / self.tg.t5) * (tg2 + a1 * tg1)

    def derivatives(self, x, v):
        n0 = self.base.n_states
        slip = float(x[1])
        tg1, tg2, tg3 = float(x[n0]), float(x[n0 + 1]), float(x[n0 + 2])
        t = self.tg
        a1 = t.t3 / t.tc
        demand = self.pref - t.kg * slip
        dev = replace(self.base, parameters=replace(self.base.parameters, tm=self.pm(x)))
        return np.concatenate([
            dev.derivatives(np.asarray(x[:n0]), v),
            [(demand - tg1) / t.ts, ((1.0 - a1) * tg1 - tg2) / t.tc, ((tg2 + a1 * tg1) * (1.0 - t.t4 / t.t5) - tg3) / t.t5],
        ])

    def injection(self, x, v):
        return self.base.injection(np.asarray(x[: self.base.n_states]), v)

    @staticmethod
    def initial(tg: PstTg, pm0: float) -> list[float]:
        a1 = tg.t3 / tg.tc
        return [pm0, (1.0 - a1) * pm0, (1.0 - tg.t4 / tg.t5) * pm0]


# ------------------------------------------------------------------- network --
def build_ybus68(scale: dict | None = None, outage: set | None = None) -> np.ndarray:
    """Ybus of the 68-bus JSON with whole-two-port branch scaling (TX4 _build_ybus convention)."""
    payload = network_payload()
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {bus: i for i, bus in enumerate(order)}
    y = np.zeros((len(order), len(order)), dtype=np.complex128)
    for e, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0.0:
            continue
        rho = (scale or {}).get(e, 1.0)
        if outage and e in outage:
            rho = 0.0
        if rho == 0.0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        series = 1.0 / complex(line["r"], line["x"])
        charging = complex(line["g"], line["b"]) / 2.0
        tap = float(line["tap"]) if float(line["tap"]) != 0.0 else 1.0
        m = tap * np.exp(1j * math.radians(float(line["phi"])))
        m2 = float(abs(m) ** 2)
        y[f, f] += rho * (series + charging) / m2
        y[t, t] += rho * (series + charging)
        y[f, t] += rho * (-series / np.conj(m))
        y[t, f] += rho * (-series / m)
    for sh in payload.get("shunts", []) or []:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return y


def branches() -> list[dict]:
    return [{"e": e, "f": int(line["bus1"]), "t": int(line["bus2"]), "r": float(line["r"]), "x": float(line["x"]),
             "tap": float(line["tap"]), "u": float(line["u"])} for e, line in enumerate(network_payload()["lines"])]


def network_with(scale: dict | None = None, outage: set | None = None):
    net = base_network()
    if not scale and not outage:
        return net
    return replace(net, ybus=build_ybus68(scale, outage))


# --------------------------------------------------------------------- cases --
@dataclass
class Case68:
    dae: object
    equilibrium: object
    gz_condition: float
    variant: str
    members: tuple
    meta: dict


DRAW_MACHINE = ("h", "xd1q1", "ka", "pss_k")
DRAW_CONV = ("pll", "current", "outer")


def _machine_params(p, variant: str, bus: int, mdraw: dict | None):
    d = pst_damping(bus) if variant in ("REAL", "NOGOV") else 0.0
    p = replace(p, d=d)
    if mdraw:
        if "h" in mdraw:
            p = replace(p, h=p.h * mdraw["h"])
        if "xd1q1" in mdraw:
            p = replace(p, xd1=p.xd1 * mdraw["xd1q1"], xq1=p.xq1 * mdraw["xd1q1"])
        if "ka" in mdraw:
            p = replace(p, gain_scale=p.gain_scale * mdraw["ka"])
        if "pss_k" in mdraw and p.pss:
            pd = dict(p.pss)
            pd["k"] = pd["k"] * mdraw["pss_k"]
            p = replace(p, pss=tuple(sorted(pd.items())))
    return p


def converter_params(g: float, cdraw: dict | None = None) -> ConverterParameters:
    c = ConverterParameters(voltage_control=True, voltage_gain=float(g), voltage_leak=LEAK)
    if cdraw:
        kw = {}
        if "pll" in cdraw:
            kw.update(kp_pll=c.kp_pll * cdraw["pll"], ki_pll=c.ki_pll * cdraw["pll"])
        if "current" in cdraw:
            kw.update(kp_i=c.kp_i * cdraw["current"], ki_i=c.ki_i * cdraw["current"])
        if "outer" in cdraw:
            kw.update(kp_p=c.kp_p * cdraw["outer"], ki_p=c.ki_p * cdraw["outer"], kp_q=c.kp_q * cdraw["outer"], ki_q=c.ki_q * cdraw["outer"])
        c = replace(c, **kw)
    return c


def build68(members, *, variant: str = "REAL", g: float = 0.0, k: float = 1.0, draw: dict | None = None, network=None, tol: float = 1e-9) -> Case68:
    """Build and equilibrate one portfolio. ``draw`` = {'machine': {...}, 'converter': {...}} multipliers."""
    assert variant in VARIANTS
    net = network or base_network()
    plan = ReplacementPlan.of({b: 1.0 for b in members})
    draw = draw or {}
    dae, x0, z0 = build_dae(plan, network=net, converter=converter_params(g, draw.get("converter")), machine_scaling={"ka": float(k)})
    pf = dae.power_flow
    slots, states, cursor = [], [], 0
    for slot in dae.slots:
        dev = slot.device
        xs = x0[slot.start: slot.stop]
        if slot.kind == "sg":
            params = _machine_params(dev.parameters, variant, slot.bus, draw.get("machine"))
            share = pf.injection(slot.bus, net.ybus) + net.loads.get(slot.bus, 0j)
            dev, xs = Ieee68Machine(bus=slot.bus, parameters=params, weight=dev.weight).initialize(pf.at(slot.bus), share)
            if variant == "REAL":
                tg = pst_tg(slot.bus)
                pm0 = float(dev.parameters.tm)
                dev = Gov68(base=dev, tg=tg, pref=pm0)
                xs = np.concatenate([xs, Gov68.initial(tg, pm0)])
        slots.append(replace(slot, device=dev, start=cursor, stop=cursor + xs.size))
        states.append(np.asarray(xs, float))
        cursor += xs.size
    new = replace(dae, slots=tuple(slots))
    new.__post_init__()
    x0n = np.concatenate(states)
    eq = solve_equilibrium(new, {}, x0n, z0, tol=tol)
    if not eq.ok:
        raise InfeasibleReplacement(f"equilibrium failed for {label(members)} ({variant}): {eq.status}")
    return Case68(dae=new, equilibrium=eq, gz_condition=float("nan"), variant=variant, members=tuple(members), meta={"g": g, "k": k})


def reduced(jac):
    return jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)


def matrices(case: Case68):
    x, z = case.equilibrium.x, case.equilibrium.z
    j1 = central_difference_jacobians(case.dae, x, z, {})
    j2 = central_difference_jacobians(case.dae, x, z, {}, scale_x=2.0, scale_z=2.0)
    a = reduced(j1)
    return a, a - reduced(j2), j1


def quotient_w(case: Case68):
    """Structural centre: span{R_x} with governors or any rotor damping; span{R_x, w} for SP33."""
    if case.variant == "SP33":
        fp = frequency_partner(case.dae)
        assert fp.w is not None, fp.reason
        return fp.w
    return None


def _participation(vals, vr, vl, j, labels):
    p = np.abs(vl[:, j].conj() * vr[:, j])
    s = p.sum()
    p = p / s if s > 0 else p
    groups = {}
    for lab, val in zip(labels, p, strict=True):
        name = lab.rsplit("_", 1)[0]
        groups[name] = groups.get(name, 0.0) + float(val)
    em = sum(v for kname, v in groups.items() if kname in ("delta", "omega"))
    return em, groups


def evaluate(case: Case68, *, modes: bool = True, participation: bool = True) -> dict:
    """TX4 direct path on the 68-bus case: transverse status/alpha, critical and EM-band modes, shapes."""
    a, d, jac = matrices(case)
    x, z = case.equilibrium.x, case.equilibrium.z
    r_x, _ = rotation_generator(case.dae, z)
    w = quotient_w(case)
    dead = np.flatnonzero(~np.any(a != 0.0, axis=1))
    keep = np.arange(a.shape[0])
    if dead.size:
        keep = np.setdiff1d(keep, dead)
        assert np.all(r_x[dead] == 0.0)
    ak, dk, rk = a[np.ix_(keep, keep)], d[np.ix_(keep, keep)], r_x[keep]
    wk = None if w is None else w[keep]
    tr = transverse_operator(ak, rk, wk)
    ev_t = np.linalg.eigvals(tr.a_perp)
    status = classify_spectrum(tr.a_perp, tr.z.T @ dk @ tr.z, SAFETY).status
    i_c = int(np.argmax(ev_t.real))
    lam_c = complex(ev_t[i_c])
    if lam_c.imag < 0:
        lam_c = lam_c.conjugate()
    others = ev_t[(np.abs(ev_t - lam_c) > 1e-9) & (np.abs(ev_t - lam_c.conjugate()) > 1e-9)]
    out = {
        "variant": case.variant, "S": label(case.members), "status": status, "alpha": float(ev_t.real.max()),
        "lam_re": lam_c.real, "lam_hz": abs(lam_c.imag) / (2 * np.pi), "rhp": int((ev_t.real > 0).sum()),
        "n_x": int(a.shape[0]), "n_dead": int(dead.size), "gap2": float(lam_c.real - others.real.max()) if others.size else float("inf"),
        "min_abs_re": float(np.abs(ev_t.real).min()), "eq_residual": float(max(case.equilibrium.norm_f, case.equilibrium.norm_g)),
        "gz_cond": float(np.linalg.cond(jac.gz)),
    }
    if modes:
        vals, vl, vr = sla.eig(a, left=True, right=True)
        phi_all = -np.linalg.solve(jac.gz, jac.gx @ vr)
        f = np.abs(vals.imag) / (2 * np.pi)
        keepm = []
        for j in range(vals.size):
            lam = vals[j]
            if abs(lam) < STRUCT_ZERO or lam.imag < -1e-12:
                continue
            if MODE_BAND[0] <= f[j] <= MODE_BAND[1] and lam.real >= MODE_RE_MIN:
                keepm.append(j)
        j_c = int(np.argmin(np.abs(vals - lam_c)))
        if j_c not in keepm:
            keepm.append(j_c)
        band = [j for j in keepm if MODE_BAND[0] <= f[j] <= MODE_BAND[1]]
        j_em = max(band, key=lambda j: vals[j].real) if band else None
        modes_out = []
        for j in keepm:
            ph = phi_all[:, j]
            n = np.linalg.norm(ph)
            ph = ph / n if n > 0 else ph
            kk = int(np.argmax(np.abs(ph)))
            ph = ph * np.exp(-1j * np.angle(ph[kk]))
            rec = {"re": float(vals[j].real), "hz": float(abs(vals[j].imag) / (2 * np.pi)), "crit": j == j_c, "em_top": j == j_em,
                   "phi": np.round(np.concatenate([ph.real, ph.imag]), 5).astype(float).tolist()}
            if participation and (j == j_c or j == j_em):
                em, groups = _participation(vals, vr, vl, j, case.dae.labels)
                rec["part_em"] = em
                rec["part_top"] = sorted(groups.items(), key=lambda kv: -kv[1])[:5]
            modes_out.append(rec)
        out["modes"] = modes_out
        out["em_top_re"] = float(vals[j_em].real) if j_em is not None else float("nan")
        out["em_top_hz"] = float(f[j_em]) if j_em is not None else float("nan")
    return out


def eval_portfolio(members, *, variant="REAL", g=0.0, k=1.0, draw=None, network=None, modes=True) -> dict:
    t0 = time.time()
    try:
        case = build68(members, variant=variant, g=g, k=k, draw=draw, network=network)
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as e:
        return {"variant": variant, "S": label(members), "status": "INFEASIBLE", "alpha": float("nan"), "error": str(e)[:300]}
    rec = evaluate(case, modes=modes)
    rec["wall_s"] = time.time() - t0
    return rec


def mac(phi1, phi2) -> float:
    a = np.asarray(phi1[: len(phi1) // 2]) + 1j * np.asarray(phi1[len(phi1) // 2:])
    b = np.asarray(phi2[: len(phi2) // 2]) + 1j * np.asarray(phi2[len(phi2) // 2:])
    num = abs(np.vdot(a, b)) ** 2
    den = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(num / den) if den > 0 else 0.0


def in_band(hz) -> bool:
    return hz is not None and np.isfinite(hz) and MODE_BAND[0] <= hz <= MODE_BAND[1]


# --------------------------------------------------------------------- store --
def task_key(task: dict) -> str:
    return hashlib.sha1(json.dumps(task, sort_keys=True, default=str).encode()).hexdigest()[:16]


class Store:
    def __init__(self, phase: str):
        self.dir = RAW / phase
        self.dir.mkdir(parents=True, exist_ok=True)

    def path(self, task):
        return self.dir / f"{task_key(task)}.json"

    def done(self, task) -> bool:
        p = self.path(task)
        if not p.exists():
            return False
        try:
            return bool(json.loads(p.read_text()).get("ok"))
        except json.JSONDecodeError:
            return False

    def put(self, task, record: dict, ok: bool = True):
        p = self.path(task)
        tmp = p.with_suffix(".tmp")
        tmp.write_text(json.dumps({"task": task, "ok": ok, "record": record}, default=str))
        tmp.replace(p)

    def all(self):
        for p in sorted(self.dir.glob("*.json")):
            try:
                yield json.loads(p.read_text())
            except json.JSONDecodeError:
                continue


def write_json(path: Path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=str), encoding="utf-8")


def log(phase: str, msg: str):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} [{phase}] {msg}"
    print(line, flush=True)
    with (LOGS / f"{phase}.log").open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
