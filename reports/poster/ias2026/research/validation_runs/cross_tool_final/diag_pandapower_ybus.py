# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Diagnostic for the Phase B static failure: where does the pandapower network differ?

Builds the pandapower net exactly as the pack does (through the compat wrapper),
extracts pandapower's internal Ybus and compares it with the project Ybus entry by
entry. Prints the branch classification made by from_ppc and the worst entries.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandapower as pp  # noqa: E402
import run_pack_pandapower_compat as compat  # noqa: E402  (applies the three fixes)

if "--tapside" in sys.argv:
    import run_pack_pandapower_tapside  # noqa: E402,F401  (rebinds compat.from_ppc)

    sys.argv.remove("--tapside")
    print("VARIANT: B-corrected (tap at canonical from bus)")
else:
    print("VARIANT: as supplied (+ compat fixes)")

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "src"))
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402

canonical = compat.canonical_ppc.canonical_to_ppc(repo)
net = compat.from_ppc(canonical, f_hz=60, validate_conversion=False)
pp.runpp(
    net,
    algorithm="nr",
    calculate_voltage_angles=True,
    init="flat",
    enforce_q_lims=False,
    tolerance_mva=1e-10,
    max_iteration=50,
    numba=False,
)

print("bus index:", list(net.bus.index)[:5], "...", len(net.bus))
print(
    "lines:",
    len(net.line),
    "trafos:",
    len(net.trafo),
    "impedances:",
    len(net.impedance),
)
lookup = net._from_ppc_lookups["branch"]
payload = json.loads((repo / "configs/ias2026/ieee39_network.json").read_text())
for k, x in enumerate(payload["lines"]):
    et = lookup.loc[k, "element_type"]
    if et != "line":
        print(
            f"  branch {k}: {int(x['bus1'])}-{int(x['bus2'])} tap={x['tap']} trans={x['trans']} "
            f"b={x['b']} -> {et} {int(lookup.loc[k, 'element'])}"
        )

ybus_pp = net._ppc["internal"]["Ybus"].toarray()
bus_lookup = net._pd2ppc_lookups["bus"]
pynet = load_network(repo / "configs/ias2026/ieee39_network.json")
# map project bus order to pandapower internal ppc order
pp_bus_for_project = []
for b in pynet.bus_idx:
    pd_idx = b if b in net.bus.index else None
    if pd_idx is None:
        raise SystemExit(f"bus {b} not found in pandapower index")
    pp_bus_for_project.append(bus_lookup[pd_idx])
pp_bus_for_project = np.asarray(pp_bus_for_project)
ybus_pp_ord = ybus_pp[np.ix_(pp_bus_for_project, pp_bus_for_project)]
diff = np.abs(ybus_pp_ord - pynet.ybus)
print("Ybus max abs diff:", diff.max())
worst = np.dstack(np.unravel_index(np.argsort(-diff.ravel())[:12], diff.shape))[0]
for i, j in worst:
    print(
        f"  ({pynet.bus_idx[i]},{pynet.bus_idx[j]}) pp={ybus_pp_ord[i, j]:.6f} "
        f"py={pynet.ybus[i, j]:.6f} |d|={diff[i, j]:.3e}"
    )
print(
    net.trafo[
        [
            "hv_bus",
            "lv_bus",
            "vn_hv_kv",
            "vn_lv_kv",
            "sn_mva",
            "vk_percent",
            "vkr_percent",
            "tap_side",
            "tap_pos",
            "tap_step_percent",
            "tap_neutral",
        ]
    ].to_string()
)
