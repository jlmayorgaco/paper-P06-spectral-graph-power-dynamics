# ruff: noqa: E501  -- sympy equation strings and parameter tables kept on one line
"""PCV06 step 1 (tx3-analysis): parameter handoff, Gate 0 vectors, internal reference.

Spec: docs/20260911_GFL_REPRODUCTION_SPEC.md (commit 3e847a4f), sections 6-7.

Writes two files so that the ANDES side never sees internal results:
- results/PCV/PCV06/PCV06_handoff_inputs.json    device parameters, policies, Gate 0
  inputs (states, terminal voltages, setpoints). Read by the ANDES side.
- results/PCV/PCV06/PCV06_internal_reference.json internal Gate 0 outputs and the 32
  internal cases (verdicts, alpha, spectra, equilibrium voltages). Read only by the
  comparison step.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from _pcv import H4, POINTS, SUBSETS4, label, out_dir, transverse_of  # noqa: E402

from _f7_common import LEAK, Theta, solve_subset  # noqa: E402
from ibr_cycles.certification.physical import physical_matrices  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import _controller_payload  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402

OUT = out_dir("PCV06")
SEED = 20260918
N_VEC = 20
CASE_POINTS = ("P4", "G_S")
GFL_FIELDS = (
    "kp_pll",
    "ki_pll",
    "tau_p",
    "kp_p",
    "ki_p",
    "kp_q",
    "ki_q",
    "kp_i",
    "ki_i",
    "xf",
    "rf",
    "kp_v",
    "ki_v",
)


def cplx(z):
    return [float(np.real(z)), float(np.imag(z))]


def transverse_spectrum(case):
    a, _, _ = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(a, r_x, frequency_partner(case.dae).w)
    return np.linalg.eigvals(tr.a_perp)


def main() -> int:
    p4 = Theta(*POINTS["P4"])
    base = solve_subset((), p4)
    h4 = solve_subset(H4, p4)
    net = base.dae.network
    payload = _controller_payload(net.config_path or None)

    # ---------------- device parameters --------------------------------------------
    sg = {}
    for bus in net.generator_buses:
        raw = net.machines[bus]
        avr = payload["avr_by_bus"][bus]
        pss = payload["pss_by_bus"][bus]
        assert float(avr.get("TB", 0.0)) == 0.0
        sg[str(bus)] = {
            "w": float(raw["Sn"]) / 100.0,
            "ra": raw["ra"],
            "xd": raw["xd"],
            "xq": raw["xq"],
            "xd1": raw["xd1"],
            "xq1": raw["xq1"],
            "Td10": raw["Td10"],
            "Tq10": raw["Tq10"],
            "M": raw["M"],
            "D": raw["D"],
            "KA_case": avr["KA"],
            "TE_case": avr["TE"],
            "KS": pss["KS"],
            "T4": pss["T4"],
            "T5": pss["T5"],
            "T6": pss["T6"],
        }
    dflt = ConverterParameters()
    gfl = {
        "defaults": {f: getattr(dflt, f) for f in GFL_FIELDS},
        "leak": LEAK,
        "w": {str(b): float(net.machines[b]["Sn"]) / 100.0 for b in H4},
    }
    # cross-check against the instantiated devices
    for slot in base.dae.slots:
        p = slot.device.parameters
        s = sg[str(slot.bus)]
        assert (
            abs(p.ka - s["KA_case"] * 1.425) < 1e-12
            and abs(p.ta - s["TE_case"] * 1.5) < 1e-12
        )
        assert abs(slot.weight - s["w"]) < 1e-12
    for slot in h4.dae.slots:
        if slot.kind == "gfl":
            assert abs(slot.weight - gfl["w"][str(slot.bus)]) < 1e-12

    # ---------------- Gate 0 vectors -----------------------------------------------
    rng = np.random.default_rng(SEED)
    x_base = base.equilibrium.x
    x_h4 = h4.equilibrium.x
    instances = [
        ("SG", slot, x_base[slot.start : slot.stop], slot.device)
        for slot in base.dae.slots
    ]
    gfl_slots = [slot for slot in h4.dae.slots if slot.kind == "gfl"]
    for g in (0.03625, 0.25):
        for slot in gfl_slots:
            dev = replace(
                slot.device, parameters=replace(slot.device.parameters, voltage_gain=g)
            )
            instances.append((f"GFL_g{g}", slot, x_h4[slot.start : slot.stop], dev))
    gate_in, gate_out = [], []
    for kind, slot, xeq, dev in instances:
        for k in range(N_VEC):
            x = xeq + 0.05 * np.maximum(1.0, np.abs(xeq)) * rng.standard_normal(
                xeq.size
            )
            v = rng.uniform(0.9, 1.1)
            a = rng.uniform(-np.pi, np.pi)
            volt = v * np.exp(1j * a)
            f = dev.derivatives(x, volt)
            s_inj = volt * np.conj(dev.injection(x, volt))
            p = dev.parameters
            setp = (
                {"pm": p.pm, "vref": p.vref}
                if kind == "SG"
                else {
                    "p_ref": p.p_ref,
                    "q_ref": p.q_ref,
                    "v_ref": p.v_ref,
                    "g": p.voltage_gain,
                }
            )
            gate_in.append(
                {
                    "kind": kind,
                    "bus": slot.bus,
                    "k": k,
                    "x": list(map(float, x)),
                    "v": float(v),
                    "a": float(a),
                    "setpoints": {n: float(val) for n, val in setp.items()},
                }
            )
            gate_out.append(
                {
                    "kind": kind,
                    "bus": slot.bus,
                    "k": k,
                    "f": list(map(float, f)),
                    "P_inj": float(s_inj.real),
                    "Q_inj": float(s_inj.imag),
                }
            )

    # ---------------- internal reference: 32 cases ---------------------------------
    truth = pd.read_csv(OUT.parent / "PCV02" / "PCV02_truth_core.csv")
    cases = []
    for pid in CASE_POINTS:
        th = Theta(*POINTS[pid])
        for s in SUBSETS4:
            case = solve_subset(s, th)
            tr = transverse_of(case)
            ev = transverse_spectrum(case)
            ref = truth[(truth.point == pid) & (truth.subset == label(s))].iloc[0]
            assert ref.status == tr["status"] and abs(ref.alpha - tr["alpha"]) < 1e-12
            volts = case.dae.voltages(case.equilibrium.z)
            cases.append(
                {
                    "point": pid,
                    "subset": label(s),
                    "status": tr["status"],
                    "alpha": tr["alpha"],
                    "rhp": tr["rhp"],
                    "crit_re": tr["crit_re"],
                    "crit_hz": tr["crit_hz"],
                    "n_x": int(case.dae.n_x),
                    "spectrum": [
                        cplx(z) for z in sorted(ev, key=lambda z: (-z.real, z.imag))
                    ],
                    "voltages": {
                        str(b): cplx(volts[case.dae.network.position(b)])
                        for b in net.bus_idx
                    },
                }
            )

    inputs = {
        "spec_commit": "3e847a4f",
        "prereg_commit": "5d0b1986",
        "system_base_mva": 100.0,
        "frequency_hz": 60.0,
        "policies": {
            pid: dict(zip("gkth", POINTS[pid], strict=True)) for pid in CASE_POINTS
        },
        "sg": sg,
        "gfl": gfl,
        "subsets": [label(s) for s in SUBSETS4],
        "gate0_inputs": gate_in,
        "gate0_seed": SEED,
    }
    reference = {"gate0_outputs": gate_out, "cases": cases}
    (OUT / "PCV06_handoff_inputs.json").write_text(
        json.dumps(inputs, indent=1), encoding="utf-8"
    )
    (OUT / "PCV06_internal_reference.json").write_text(
        json.dumps(reference, indent=1), encoding="utf-8"
    )
    print(
        "instances", len(instances), "gate0 vectors", len(gate_in), "cases", len(cases)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
