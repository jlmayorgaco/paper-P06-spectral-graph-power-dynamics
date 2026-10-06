"""Frozen source and generator-component share invariants for Experiment N."""
from __future__ import annotations

import csv
import hashlib
import math
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


class FrozenContract(unittest.TestCase):
    def test_environment_and_n01_are_frozen(self) -> None:
        provenance = tomllib.loads((OUT / "SOFTWARE_PROVENANCE.toml").read_text())
        self.assertTrue(provenance["same_environment_as_expM"])
        for name in ("Project", "Manifest"):
            self.assertEqual(hashlib.sha256((ROOT / f"{name}.toml").read_bytes()).hexdigest(),
                             provenance[f"{name.lower()}_sha256"])
        path = OUT / "TABLE_N01_original_operating_point.csv"
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         (OUT / (path.name + ".sha256")).read_text().strip())

    def test_component_dispatch_and_separate_zip_load(self) -> None:
        with (OUT / "TABLE_N01_original_operating_point.csv").open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual([int(r["bus"]) for r in rows], list(range(30, 40)))
        total = sum(float(r["P_gen_MW"]) for r in rows)
        self.assertTrue(math.isclose(total, 5402.761089978776, abs_tol=1e-8))
        for r in rows:
            self.assertLess(abs(float(r["P_gen_MW"]) + float(r["P_load_MW"])
                                - float(r["P_net_MW"])), 1e-8)
            self.assertLess(abs(float(r["Q_gen_Mvar"]) + float(r["Q_load_Mvar"])
                                - float(r["Q_net_Mvar"])), 1e-8)
        self.assertNotEqual(float(rows[1]["P_load_MW"]), 0.0)
        self.assertNotEqual(float(rows[9]["P_load_MW"]), 0.0)

    def test_rho_contract_conserves_component_power(self) -> None:
        with (OUT / "TABLE_N02_rho_semantics.csv").open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 60)
        for r in rows:
            rho = float(r["rho"])
            self.assertTrue(math.isclose(float(r["SG_target_P_MW"]) +
                                          float(r["GFL_target_P_MW"]),
                                          float(r["total_target_P_MW"]), abs_tol=1e-9))
            self.assertTrue(math.isclose(float(r["SG_target_Q_Mvar"]) +
                                          float(r["GFL_target_Q_Mvar"]),
                                          float(r["total_target_Q_Mvar"]), abs_tol=1e-9))
            self.assertTrue(math.isclose(float(r["GFL_port_scale"]), rho, abs_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
