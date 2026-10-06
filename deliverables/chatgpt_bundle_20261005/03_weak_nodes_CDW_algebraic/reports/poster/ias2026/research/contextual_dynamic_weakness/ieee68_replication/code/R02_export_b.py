# ruff: noqa: E501
"""R2 Model B handoff (tx3-analysis): everything the ANDES side needs, taken from the instantiated internal
devices (no internal result other than the power flow). Writes results/b68/B68_handoff.json.

- network rows (buses, lines, PV, slack, loads) from configs/ieee68/ieee68_network.json;
- per-bus SG68 parameters (exact internal values, incl. the DC4B Aex/Bex and the KI/KD mapping of
  Init_MultiMachine.m) and the PST governor / damping of the REAL variant;
- converter targets: base power-flow P, Q of every candidate bus and the rating |S|/0.8 (x 100 MVA);
- the frozen WECC library configuration (TX3-GFL-0.1) as embedded in the hardening handoff (sha f07a6a40).
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json

B68 = R.RESULTS / "b68"
B68.mkdir(parents=True, exist_ok=True)


def main():
    R.env_threads()
    case = R.build68((), variant="REAL", g=0.0, k=1.0)
    net = R.base_network()
    pf = case.dae.power_flow
    sg = {}
    for s in case.dae.slots:
        dev = s.device
        p = dev.parameters
        rec = {"cls": {"DC4B": "SG68D", "ST1A": "SG68S", "MANUAL": "SG68M"}[p.exciter], "w": float(dev.weight),
               "xl": p.xl, "ra": p.ra, "xd": p.xd, "xd1": p.xd1, "xd2": p.xd2, "Td10": p.td10, "Td20": p.td20,
               "xq": p.xq, "xq1": p.xq1, "xq2": p.xq2, "Tq10": p.tq10, "Tq20": p.tq20, "H": p.h, "D": p.d, "Tcd": p.tc,
               "KG": dev.tg.kg, "TSG": dev.tg.ts, "TCG": dev.tg.tc, "T3G": dev.tg.t3, "T4G": dev.tg.t4, "T5G": dev.tg.t5}
        e = p.e
        if p.exciter == "DC4B":
            rec.update({"TR": e["tr"], "KA": e["ka"], "TA": e["ta"], "KE": e["ke"], "TE": e["te"], "AEX": e["aex"], "BEX": e["bex"],
                        "KF": e["kf"], "TF": e["tf"], "KP": e["kp"], "KI": e["ki"], "KD": e["kd"], "TD": e["td"]})
        elif p.exciter == "ST1A":
            rec.update({"TR": e["tr"], "KA": e["ka"]})
        if p.pss:
            q = p.p
            rec.update({"K": q["k"], "TW": q["tw"], "T11": q["t11"], "T12": q["t12"], "T21": q["t21"], "T22": q["t22"], "T31": q["t31"], "T32": q["t32"]})
        sg[str(s.bus)] = rec
    targets = {}
    for b in R.PHYSICAL:
        gen = pf.injection(b, net.ybus) + net.loads.get(b, 0j)
        targets[str(b)] = {"p": float(gen.real), "q": float(gen.imag), "rating_mva": float(abs(gen) / net.converter_loading * 100.0)}
    payload = R.network_payload()
    hard = json.loads((R.CDW / "results" / "hardening" / "alt" / "H17_handoff.json").read_text(encoding="utf-8"))
    design = json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())
    pv_v = {int(b): abs(pf.at(b)) for b in net.generator_buses}
    handoff = {"prereg_commit": "6f188531", "network": {k: payload[k] for k in ("buses", "lines", "pv", "slack", "loads")},
               "pf_voltage": {str(b): v for b, v in pv_v.items()},
               "sg": sg, "targets": targets, "wecc_config": hard["wecc_config"], "wecc_config_sha256": hard["wecc_config_sha256"],
               "policies_B": design["policies_B"], "draws_B": design["draws_B"], "V68": list(R.V68), "r11_branches": design["r11_branches"],
               "eligible_branches": design["eligible_branches"], "gammas": design["gammas"], "ybus_internal_re": R.build_ybus68().real.tolist(),
               "ybus_internal_im": R.build_ybus68().imag.tolist()}
    (B68 / "B68_handoff.json").write_text(json.dumps(handoff, indent=1, default=str), encoding="utf-8")
    print("sg", len(sg), "targets", len(targets), "sha", handoff["wecc_config_sha256"][:12])
    assert handoff["wecc_config_sha256"].startswith("f07a6a40")


if __name__ == "__main__":
    main()
