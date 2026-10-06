"""Materialize a one-gain uncertainty case without changing the parent design."""
from __future__ import annotations

import sys
import tomllib
from pathlib import Path

here = Path(__file__).resolve().parent
base_id, output_id, kind, bus, factor = sys.argv[1:]
bus = int(bus)
factor = float(factor)
assert kind in ("Kp", "Ki") and 30 <= bus <= 39 and factor > 0
base = tomllib.loads((here / "designs" / f"{base_id}.toml").read_text())
kp, ki = list(base["Kp"]), list(base["Ki"])
(kp if kind == "Kp" else ki)[bus - 30] *= factor
out = here / "designs" / f"{output_id}.toml"
assert not out.exists(), f"refusing to overwrite {out}"
out.write_text(
    "rho = " + repr(base["rho"]) + "\n"
    + "Kp = " + repr(kp) + "\n"
    + "Ki = " + repr(ki) + "\n"
    + f'uncertainty_parent = "{base_id}"\n'
    + f'gain_kind = "{kind}"\nperturbed_bus = {bus}\nfactor = {factor}\n',
    encoding="utf-8",
)
print(out)
