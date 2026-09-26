"""Export frozen IEEE-39 JSON rows as dependency-light Julia constants."""

from __future__ import annotations

import json
from pathlib import Path

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
CONFIG = ROOT / "raw" / "true_same_model" / "canonical_source" / "configs" / "ieee39_network.json"
OUTPUT = ROOT / "code" / "julia" / "frozen_ieee39_data.jl"


def f(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return repr(float(value)) if isinstance(value, (int, float)) else repr(value)


def rows(items: list[dict], fields: tuple[str, ...]) -> str:
    return "[" + ",\n".join("(" + ", ".join(f(item[field]) for field in fields) + ")" for item in items) + "]"


payload = json.loads(CONFIG.read_text(encoding="utf-8"))
machine_fields = ("bus", "Sn", "D", "M", "ra", "xd", "xq", "xd1", "xq1", "Td10", "Tq10")
avr_fields = ("TA", "KA", "TE")
pss_fields = ("KS", "T5", "T6", "T4")
pv_fields = ("bus", "p0", "q0", "v0", "pmax", "pmin", "qmax", "qmin")
line_fields = ("bus1", "bus2", "r", "x", "g", "b", "tap", "phi", "u")
load_fields = ("bus", "p0", "q0")
shunt_fields = ("bus", "g", "b")

machines = payload["machines"]
avrs = payload["avr"]
psss = payload["pss"]
for machine, avr, pss in zip(machines, avrs, psss, strict=True):
    machine["bus"] = int(machine["bus"])
    avr["TA"] = avr.get("TA", 0.0)
    avr["KA"] = avr.get("KA", 0.0)
    avr["TE"] = avr.get("TE", 0.0)
    pss["KS"] = pss.get("KS", 0.0)
    pss["T5"] = pss.get("T5", 1.0)
    pss["T6"] = pss.get("T6", 1.0)
    pss["T4"] = pss.get("T4", 1.0)

text = """# Frozen IEEE-39 data exported from canonical_source/configs/ieee39_network.json.
# This file is generated from the frozen JSON; do not hand-edit numeric rows.
const SYSTEM_BASE_MVA = 100.0
const FBASE_HZ = 60.0
const OMEGA_B = 2.0 * pi * FBASE_HZ
const BUS_IDS = Int[%s]
const BUS_BASE_KV = Float64[%s]
const MACHINE_ROWS = %s
const AVR_ROWS = %s
const PSS_ROWS = %s
const PV_ROWS = %s
const LOAD_ROWS = %s
const LINE_ROWS = %s
const SHUNT_ROWS = %s
const SLACK_ROW = (%s, %s, %s, %s)
""" % (
    ", ".join(str(int(row["idx"])) for row in payload["buses"]),
    ", ".join(f(row["Vn"]) for row in payload["buses"]),
    rows(machines, machine_fields),
    rows(avrs, avr_fields),
    rows(psss, pss_fields),
    rows(payload["pv"], pv_fields),
    rows(payload["loads"], load_fields),
    rows(payload["lines"], line_fields),
    rows(payload["shunts"], shunt_fields),
    f(payload["slack"][0]["bus"]), f(payload["slack"][0]["v0"]),
    f(payload["slack"][0]["a0"]), f(payload["slack"][0]["Sn"]),
)
OUTPUT.write_text(text, encoding="utf-8")
print(OUTPUT)
