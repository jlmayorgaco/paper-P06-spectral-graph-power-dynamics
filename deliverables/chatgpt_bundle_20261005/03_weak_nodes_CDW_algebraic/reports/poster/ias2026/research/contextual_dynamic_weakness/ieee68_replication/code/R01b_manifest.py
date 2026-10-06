# ruff: noqa: E501
"""R1: results/CDW68_MODEL_MANIFEST.csv - every model component with its value, source and hash."""

from __future__ import annotations

import _r68 as R  # noqa: I001

import hashlib
import json

import pandas as pd


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    net = R.network_payload()
    audit = json.loads((R.RESULTS / "CDW68_R01_audit.json").read_text())
    real = next(v for v in audit["variants"] if v["variant"] == "REAL")
    rows = [
        ("network", "IEEE 68-bus NETS-NYPS, 68 buses, 83 branch rows", net["source"]["document"], sha(R.NETWORK)),
        ("buses", len(net["buses"]), "ieee68_network.json", ""),
        ("branches", len(net["lines"]), "ieee68_network.json (header says 86)", ""),
        ("transformer_rows", sum(1 for ln in net["lines"] if float(ln["tap"]) not in (0.0, 1.0)), "ieee68_network.json tap column", ""),
        ("base_mva", 100, "Singh & Pal v3.3", ""),
        ("frequency_hz", net["model"]["frequency_hz"], "Singh & Pal v3.3", ""),
        ("machines", len(net["machines"]), "sub-transient (4 rotor coils + dummy coil), eqs (1)-(6)", sha(R.RESEARCH / "src/ibr_cycles/models/ieee68_devices.py")),
        ("machine_ratings", "PST mac_con col 3: 300,800,800,800,700,900,800,800,1000,1200,1600,1900,12000,10000,10000,11000 MVA", "PST data16m.m (GridSTAGE)", sha(R.INPUTS / "source_excerpts/gridstage_pst_data16m.m")),
        ("damping_REAL", "PST d_o: 0 on G1-G12; 4.0782, 3, 3, 4.45 pu (rating base) on G13-G16", "PST data16m.m", ""),
        ("damping_SP33", "D = 0 all machines", "Singh & Pal v3.3", ""),
        ("governors_REAL", "PST tg model 1 on every surviving machine: 1/R=25 (rating base), Ts 0.1, Tc 0.5, T3 0, T4 1.25, T5 5.0 s", "PST data16m.m tg_con (commented) + pstess/tg.m", sha(R.INPUTS / "source_excerpts/pstess_tg.m")),
        ("governors_SP33", "none (report ignores governors)", "Singh & Pal v3.3", ""),
        ("avr", "DC4B on G1-G8, G10-G12; ST1A on G9; manual on G13-G16", "Singh & Pal v3.3 eqs (7)-(8)", ""),
        ("pss", "speed washout + 3 lead-lags on G1-G12", "Singh & Pal v3.3 eq (9)", ""),
        ("loads", "constant impedance from the power-flow P/Q at the solved voltage", "Gate 3 convention", ""),
        ("slack", "G16 (bus 16)", "Singh & Pal v3.3", ""),
        ("equilibrium_residual_REAL_base", real["eq_residual"], "code/R01_audit.py", ""),
        ("structural_zeros_REAL_base", real["n_zero_full_A"], "code/R01_audit.py", ""),
        ("power_flow_reproduction", json.dumps(audit["gate3_power_flow"]), "frozen G3 eligibility()", ""),
        ("converter_A", "TX4 GFL, ConverterParameters defaults, leaky Q/V (leak 0.05), rating |S|/0.8", "TX4 freeze 69f200df", sha(R.RESEARCH / "src/ibr_cycles/models/ieee39_devices.py")),
        ("converter_B", "WECC library chain TX3-GFL-0.1 (PLL2, BusFreq, REGCP1, REECB1, REPCA1), constant Q", "H17_andes_alt.py / TX3 freeze", ""),
        ("primary_frequency_inputs", "variants REAL / NOGOV / SP33", "inputs/pst_primary_frequency_v1.json", sha(R.PRIMARY)),
    ]
    pd.DataFrame(rows, columns=["component", "value", "source", "sha256"]).to_csv(R.RESULTS / "CDW68_MODEL_MANIFEST.csv", index=False)
    print(len(rows), "rows")


if __name__ == "__main__":
    main()
