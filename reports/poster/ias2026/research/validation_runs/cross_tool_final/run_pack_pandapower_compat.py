"""Run the supplied pack script ``run_pandapower_parity.py`` byte-for-byte unchanged.

Two defects stop the supplied pack from running against its own recommended stack
(pandapower 3.4.x and the frozen canonical JSON). This wrapper fixes exactly those
two, in memory, and changes nothing else (no data, no tolerances, no solver options):

1. ``from pandapower.converter import from_ppc`` fails in pandapower 3.4.0, where
   ``pandapower/converter/__init__.py`` is empty. The wrapper binds the same function,
   ``pandapower.converter.pypower.from_ppc``, under the name the script imports.
2. ``canonical_ppc.canonical_to_ppc`` reads the PV-generator keys ``p``, ``v`` and
   ``sn``. The canonical JSON stores them as ``p0``, ``v0`` and ``Sn``; the project
   loader ``ibr_cycles.models.ieee39_network.load_network`` maps them one-to-one
   (lines 126-137). The wrapper adds the three aliases with identical values. No
   value is modified.
3. Upstream pandapower 3.4.0 bug, ``converter/pypower/from_ppc.py`` line 303: for an
   "impedance" branch (tap 1.0 between different voltage levels) with RATE_A = 0 it
   indexes the transformer array ``sn`` instead of ``sn_mva`` and raises IndexError.
   The pack writes RATE_A = 0 for every branch. pandapower's own intended fallback
   for a zero rating is MAX_VAL = 99999 (lines 199, 215, 250, 303), so the wrapper
   writes that same value into RATE_A for zero-rated branches before conversion.
   RATE_A is a thermal rating and the per-unit base of the converted element; the
   electrical parameters are rescaled consistently, so the admittances are unchanged.
   The supplementary check verifies this by comparing pandapower's internal Ybus with
   the project Ybus.

Usage: python run_pack_pandapower_compat.py <research_root>
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

import pandapower.converter as converter
from pandapower.converter.pypower import from_ppc

HERE = Path(__file__).resolve().parent
PACK = HERE / "pack"

converter.from_ppc = from_ppc

sys.path.insert(0, str(PACK))
import canonical_ppc  # noqa: E402

_original_load_payload = canonical_ppc.load_payload


def _load_payload_with_aliases(repo: Path) -> dict:
    payload = _original_load_payload(repo)
    for row in payload["pv"]:
        row.setdefault("p", row["p0"])
        row.setdefault("v", row["v0"])
        row.setdefault("sn", row["Sn"])
    return payload


canonical_ppc.load_payload = _load_payload_with_aliases

PANDAPOWER_MAX_VAL = 99999.0  # from_ppc.py line 199
_original_canonical_to_ppc = canonical_ppc.canonical_to_ppc


def _canonical_to_ppc_rated(repo: Path) -> dict:
    ppc = _original_canonical_to_ppc(repo)
    rate_a = ppc["branch"][:, canonical_ppc.RATE_A]
    rate_a[rate_a == 0.0] = PANDAPOWER_MAX_VAL
    return ppc


canonical_ppc.canonical_to_ppc = _canonical_to_ppc_rated

if __name__ == "__main__":
    runpy.run_path(str(PACK / "run_pandapower_parity.py"), run_name="__main__")
