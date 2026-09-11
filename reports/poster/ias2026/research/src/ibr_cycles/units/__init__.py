"""Canonical physical quantities of a replacement portfolio (unit data contract).

See quantities.py. The legacy ``ReplacementCase.replaced_mw`` is left untouched
for the reproducibility of frozen results; it is an apparent-power rating in MVA
and must not be read as megawatts.
"""

from .quantities import (
    PMAX_DOCUMENTED,
    ReplacementQuantities,
    measure,
    per_machine_table,
)

__all__ = ["PMAX_DOCUMENTED", "ReplacementQuantities", "measure", "per_machine_table"]
