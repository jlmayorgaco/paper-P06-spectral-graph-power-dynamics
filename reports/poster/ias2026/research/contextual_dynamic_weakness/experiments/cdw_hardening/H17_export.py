# ruff: noqa: E501
"""H17 step 1 (tx3-analysis): handoff for the ALT-WECC model (prereg H17/H18).

Writes results/hardening/alt/H17_handoff.json (read by the ANDES side; no internal results) and
results/hardening/alt/H17_internal_reference.json (internal base alpha/status for gate Q1; read
only by the comparison step). SG2AX parameters are taken from the instantiated internal devices
at each policy (exact per-bus K_A and T_A, including the h heterogeneity)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import hashlib
import json

import yaml

import _cdw as C
from ibr_cycles.models.ieee39_case import _controller_payload

import H14_topology as H14
import H09_corridors as H9
import E07_topology as E7

OUT = HI.RESULTS / "alt"
OUT.mkdir(parents=True, exist_ok=True)
CFG = C.PROJECT.parents[4] / "configs" / "tx3_gfl_common.yaml"
LINES12 = (3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44)  # cross-tool holdout (seed 20260911)
CORRIDORS = ("TXother", "TXall", "K2_01")


def sg_params(theta):
    case = C.solve((), theta)
    net = case.dae.network
    payload = _controller_payload(net.config_path or None)
    out = {}
    for slot in case.dae.slots:
        p = slot.device.parameters
        raw = net.machines[slot.bus]
        pss = payload["pss_by_bus"][slot.bus]
        assert abs(p.pss_gain - pss["KS"]) < 1e-12
        out[str(slot.bus)] = {"w": float(raw["Sn"]) / 100.0, "ra": p.ra, "xd": p.xd, "xq": p.xq, "xd1": p.xd1, "xq1": p.xq1,
                              "Td10": p.td10, "Tq10": p.tq10, "M": p.m, "D": p.d, "KA": p.ka, "TE": p.ta,
                              "KS": pss["KS"], "T4": pss["T4"], "T5": pss["T5"], "T6": pss["T6"]}
        # internal mapping (ieee39_case.py): washout = T5, wash_lag = T6, lag = T4
        assert abs(p.pss_washout - pss["T5"]) < 1e-12 and abs(p.pss_wash_lag - pss["T6"]) < 1e-12 and abs(p.pss_lag - pss["T4"]) < 1e-12
    return out


def networks():
    nets = {"NOMINAL": {"scale": {}, "outage": []}}
    for e in LINES12:
        for g in (1.5, 1.001, 0.999):
            nets[f"L{e}_g{g}"] = {"scale": {str(e): g}, "outage": []}
    cor = H9.corridor_defs()
    for c in CORRIDORS:
        nets[f"C_{c}"] = {"scale": {str(k): v for k, v in H9.scale_of(cor[c], "L1_0.50").items()}, "outage": []}
    acts = {a["action"]: a for a in E7.actions()}
    for a in H14.H18_ACTIONS:
        nets[f"T_{a}"] = {"scale": acts[a]["scale"], "outage": acts[a]["outage"]}
    return nets


def main():
    pols = {p: list(HI.theta_of(p)) for p in H14.H18_POLICIES}
    cfg_bytes = CFG.read_bytes()
    qual = {"D01": list(C.policy("D01"))}  # Q3 at P4 (prereg H17)
    sg = {p: sg_params(tuple(t)) for p, t in {**pols, **qual}.items()}
    handoff = {"prereg_commit": HI.PREREG_COMMIT, "policies": pols, "qual_policies": qual, "sg": sg,
               "wecc_config": yaml.safe_load(cfg_bytes.decode("utf-8")), "wecc_config_sha256": hashlib.sha256(cfg_bytes).hexdigest(),
               "subsets": [C.label(s) for s in C.subsets(C.V4)], "networks": networks(), "lines12": list(LINES12),
               "corridors": list(CORRIDORS), "actions": list(H14.H18_ACTIONS)}
    (OUT / "H17_handoff.json").write_text(json.dumps(handoff, indent=1, default=str), encoding="utf-8")
    ref = {}
    for p, t in pols.items():
        r = C.eval_portfolio((), tuple(t), modes=False)
        ref[p] = {"alpha": r["alpha"], "status": r["status"], "lam_hz": r["lam_hz"]}
    (OUT / "H17_internal_reference.json").write_text(json.dumps(ref, indent=1), encoding="utf-8")
    print("policies", list(pols), "networks", len(handoff["networks"]), "sha", handoff["wecc_config_sha256"][:12])
    assert handoff["wecc_config_sha256"].startswith("f07a6a40")


if __name__ == "__main__":
    main()
