# ruff: noqa: E501  -- sympy equation strings and parameter tables kept on one line
"""PCV06 step 2 (xtool-andes-gfl ONLY): ANDES Gate 0 evaluation and the 32 cases.

Spec: docs/20260911_GFL_REPRODUCTION_SPEC.md (commit 3e847a4f). Reads only
PCV06_handoff_inputs.json (never the internal reference). Launch through
launch_andes.ps1 / with USERPROFILE set to .venv/xtool-andes-gfl/home.

    python run_andes.py gate0   -> results/PCV/PCV06/PCV06_andes_gate0.json
    python run_andes.py cases   -> results/PCV/PCV06/PCV06_andes_cases.json
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
assert importlib.util.find_spec("ibr_cycles") is None, (
    "internal package must not be importable"
)

import andes  # noqa: E402

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
VENV = REPO / ".venv" / "xtool-andes-gfl"
HOME = VENV / "home"
PYCODE = HOME / "pycode"
WORK = VENV / "work"
OUT = RESEARCH / "results" / "PCV" / "PCV06"
KEEP = ("Bus", "PQ", "PV", "Slack", "Shunt", "Line", "Area")
SG_STATES = ("delta", "omega", "eq1", "ed1", "efd", "pw", "pl")
GFL_STATES = (
    "theta",
    "x_pll",
    "p_f",
    "q_f",
    "x_p",
    "x_q",
    "i_d",
    "i_q",
    "x_id",
    "x_iq",
    "x_v",
)
BAND = (0.3, 1.5)
INPUTS: dict = {}
STATIC: Path = WORK / "ieee39_static.xlsx"


def static_case() -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    target = WORK / "ieee39_static.xlsx"
    src = Path(andes.__file__).parent / "cases" / "ieee39" / "ieee39_full.xlsx"
    sheets = pd.read_excel(src, sheet_name=None)
    with pd.ExcelWriter(target) as writer:
        for name in KEEP:
            if name in sheets and len(sheets[name]):
                sheets[name].to_excel(writer, sheet_name=name, index=False)
    return target


def build(policy: dict, members: tuple[int, ...]):
    andes.config_logger(stream_level=50)
    ss = andes.load(
        str(STATIC),
        setup=False,
        no_output=True,
        default_config=True,
        pycode_path=str(PYCODE),
    )
    gen_of = {int(b): i for b, i in zip(ss.PV.bus.v, ss.PV.idx.v, strict=True)}
    gen_of.update(
        {int(b): i for b, i in zip(ss.Slack.bus.v, ss.Slack.idx.v, strict=True)}
    )
    for bus in sorted(gen_of):
        if bus in members:
            gd = INPUTS["gfl"]
            ss.add(
                "GFL11",
                {
                    "idx": f"GFL_{bus}",
                    "name": f"GFL_{bus}",
                    "u": 1,
                    "bus": bus,
                    "gen": gen_of[bus],
                    "w": gd["w"][str(bus)],
                    "g": policy["g"],
                    "leak": gd["leak"],
                    **gd["defaults"],
                },
            )
        else:
            p = INPUTS["sg"][str(bus)]
            ss.add(
                "SG2AX",
                {
                    "idx": f"SG_{bus}",
                    "name": f"SG_{bus}",
                    "u": 1,
                    "bus": bus,
                    "gen": gen_of[bus],
                    "w": p["w"],
                    "ra": p["ra"],
                    "xd": p["xd"],
                    "xq": p["xq"],
                    "xd1": p["xd1"],
                    "xq1": p["xq1"],
                    "Td10": p["Td10"],
                    "Tq10": p["Tq10"],
                    "M": p["M"],
                    "D": p["D"],
                    "KA": p["KA_case"] * policy["k"],
                    "TE": p["TE_case"] * policy["t"],
                    "KS": p["KS"],
                    "T4": p["T4"],
                    "T5": p["T5"],
                    "T6": p["T6"],
                },
            )
    ss.setup()
    ss.PQ.config.pq2z = 0
    ss.PQ.config.p2p, ss.PQ.config.q2q = 1.0, 1.0
    ss.PQ.config.p2z, ss.PQ.config.q2z = 0.0, 0.0
    ss.PQ.config.p2i, ss.PQ.config.q2i = 0.0, 0.0
    ss.PV.config.pv2pq = 0
    ss.PFlow.config.tol = 1e-12
    ss.TDS.config.criteria = 0
    ss.PFlow.run()
    assert ss.PFlow.converged
    ss.TDS.init()
    return ss


def residual(ss) -> float:
    ss.TDS.fg_update(ss.exist.tds)
    return float(max(np.abs(ss.dae.f).max(), np.abs(ss.dae.g).max()))


def evaluate(model, j, names, x, v, a, setpoints):
    inp = {
        k: np.array(val, dtype=float, copy=True)
        for k, val in model.get_inputs(refresh=True).items()
    }
    for name, value in zip(names, x, strict=True):
        inp[name][j] = value
    inp["v"][j], inp["a"][j] = v, a
    for name, value in setpoints.items():
        inp[name][j] = value
    n = model.n
    f_ret = model.calls.f(*[inp[k] for k in model.calls.f_args])
    g_ret = model.calls.g(*[inp[k] for k in model.calls.g_args])
    f_names = list(model.cache.states_and_ext)
    g_names = list(model.cache.algebs_and_ext)
    f_val = {
        nm: np.broadcast_to(np.asarray(r, dtype=float), (n,))[j]
        for nm, r in zip(f_names, f_ret, strict=True)
    }
    g_val = {
        nm: np.broadcast_to(np.asarray(r, dtype=float), (n,))[j]
        for nm, r in zip(g_names, g_ret, strict=True)
    }
    return [float(f_val[nm]) for nm in names], -float(g_val["a"]), -float(g_val["v"])


def gate0() -> dict:
    p4 = INPUTS["policies"]["P4"]
    systems = {
        "SG": build(p4, ()),
        "GFL": build(p4, tuple(int(b) for b in INPUTS["gfl"]["w"])),
    }
    rows, init_setpoints = [], {}
    for rec in INPUTS["gate0_inputs"]:
        kind = "SG" if rec["kind"] == "SG" else "GFL"
        ss = systems[kind]
        model = ss.SG2AX if kind == "SG" else ss.GFL11
        j = list(model.idx.v).index(f"{kind}_{rec['bus']}")
        names = SG_STATES if kind == "SG" else GFL_STATES
        f, p, q = evaluate(
            model, j, names, rec["x"], rec["v"], rec["a"], rec["setpoints"]
        )
        rows.append(
            {
                "kind": rec["kind"],
                "bus": rec["bus"],
                "k": rec["k"],
                "f": f,
                "P_inj": p,
                "Q_inj": q,
            }
        )
        key = f"{kind}_{rec['bus']}"
        if key not in init_setpoints:
            svc = ("pm", "vref") if kind == "SG" else ("p_ref", "q_ref", "v_ref")
            init_setpoints[key] = {s: float(getattr(model, s).v[j]) for s in svc}
    return {
        "rows": rows,
        "andes_initialized_setpoints": init_setpoints,
        "residual_base": residual(systems["SG"]),
        "residual_h4": residual(systems["GFL"]),
    }


def transverse(mu):
    order = np.argsort(np.abs(mu), kind="stable")
    small = np.abs(mu[order[:2]])
    rest = mu[order[2:]]
    ok = bool(small.max() < 1e-3 and np.abs(rest).min() >= 1e-2)
    return rest, small, ok


def case_metrics(ev):
    f = ev.imag / (2 * np.pi)
    band = ev[(f >= BAND[0]) & (f <= BAND[1])]
    crit = complex(band[np.argmax(band.real)]) if band.size else complex("nan")
    rhp = int((ev.real > 0).sum())
    resolved = float(np.abs(ev.real).min()) >= 2e-3
    verdict = (
        ("UNSTABLE" if rhp else "STABLE") if resolved else "BOUNDARY_OR_UNRESOLVED"
    )
    return {
        "alpha": float(ev.real.max()),
        "rhp": rhp,
        "crit_re": crit.real,
        "crit_hz": crit.imag / (2 * np.pi),
        "min_abs_re": float(np.abs(ev.real).min()),
        "verdict": verdict,
    }


def hypergraph(verdicts: dict) -> dict:
    if verdicts[()] != "STABLE":
        return {"H": f"BASE_{verdicts[()]}", "kappa": None}
    unsafe = [s for s, v in verdicts.items() if v == "UNSTABLE"]
    minimal = sorted(
        (s for s in unsafe if not any(set(r) < set(s) for r in unsafe)),
        key=lambda s: (len(s), s),
    )
    lab = ["+".join(map(str, e)) for e in minimal]
    return {
        "H": "|".join(lab) or "EMPTY",
        "kappa": min((len(e) for e in minimal), default=None),
        "exact": not any(v == "BOUNDARY_OR_UNRESOLVED" for v in verdicts.values()),
    }


def cases() -> dict:
    out = {"cases": [], "hypergraphs": {}}
    for pid, policy in INPUTS["policies"].items():
        verdicts = {}
        for lab in INPUTS["subsets"]:
            members = () if lab == "BASE" else tuple(int(b) for b in lab.split("+"))
            ss = build(policy, members)
            res = residual(ss)
            ss.EIG.run()
            mu = np.asarray(ss.EIG.mu)
            ev, small, pair_ok = transverse(mu)
            met = case_metrics(ev)
            verdicts[members] = met["verdict"]
            volts = {
                str(int(b)): [float(v * np.cos(a)), float(v * np.sin(a))]
                for b, v, a in zip(ss.Bus.idx.v, ss.Bus.v.v, ss.Bus.a.v, strict=True)
            }
            out["cases"].append(
                {
                    "point": pid,
                    "subset": lab,
                    **met,
                    "init_residual": res,
                    "n_eig": int(mu.size),
                    "n_x": int(ss.dae.n),
                    "structural_pair": [float(s) for s in small],
                    "structural_pair_ok": pair_ok,
                    "spectrum": [
                        [float(z.real), float(z.imag)]
                        for z in sorted(ev, key=lambda z: (-z.real, z.imag))
                    ],
                    "voltages": volts,
                }
            )
            print(
                f"  {pid:4s} {lab:12s} {met['verdict']:9s} a={met['alpha']:+.6f} f={met['crit_hz']:.5f} "
                f"rhp={met['rhp']} res={res:.1e} pair={pair_ok}",
                flush=True,
            )
        out["hypergraphs"][pid] = hypergraph(verdicts)
    return out


if __name__ == "__main__":
    assert Path(os.path.expanduser("~")).resolve() == HOME.resolve(), (
        "USERPROFILE must point to the private home"
    )
    assert Path(andes.__file__).resolve().is_relative_to((VENV / "src").resolve())
    INPUTS.update(
        json.loads((OUT / "PCV06_handoff_inputs.json").read_text(encoding="utf-8"))
    )
    STATIC = static_case()
    mode = sys.argv[1]
    result = gate0() if mode == "gate0" else cases()
    result["andes_version"] = andes.__version__
    result["andes_file"] = str(Path(andes.__file__).relative_to(REPO)).replace(
        "\\", "/"
    )
    (OUT / f"PCV06_andes_{mode}.json").write_text(
        json.dumps(result, indent=1), encoding="utf-8"
    )
    print("wrote", f"PCV06_andes_{mode}.json")
