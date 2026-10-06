"""Integrity gates for the source-of-truth inputs used in Experiment M."""
from __future__ import annotations

import csv
import hashlib
import math
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports" / "experiment_M"


class AuditIntegrity(unittest.TestCase):
    def test_provenance_matches_current_inputs(self) -> None:
        provenance = tomllib.loads((REPORT / "SOFTWARE_PROVENANCE.toml").read_text(encoding="utf-8"))
        for name in ("Project.toml", "Manifest.toml"):
            expected = provenance[name.split(".")[0].lower() + "_sha256"]
            self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), expected)
        for item in provenance["ieee39_inputs"].values():
            self.assertEqual(hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest(), item["sha256"])

    def test_component_bases_are_uniform(self) -> None:
        with (REPORT / "tables" / "TABLE_M01_component_bases.csv").open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertGreaterEqual(len(rows), 340)
        self.assertTrue(all(abs(float(r["Sbase"]) - 100) < 1e-12 for r in rows))
        self.assertTrue(all(abs(float(r["fbase_equivalent"]) - 60) < 1e-10 for r in rows))
        failures = [r for r in rows if r["pass"].lower() != "true"]
        self.assertTrue(any(r["case"] == "ExpK_nominal" and r["component"] == "bus30_GFL"
                            and float(r["Vbase"]) == 1.0 for r in failures))
        self.assertTrue(any(r["case"] == "GFL_template" for r in rows))
        self.assertTrue(any(r["case"] == "mixed_bus38" for r in rows))

    def test_frozen_case_definitions_cover_every_generator(self) -> None:
        with (REPORT / "tables" / "TABLE_M05_case_definitions.csv").open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 30)
        for case in ("all_SG", "ExpG_candidate", "ExpK_nominal"):
            self.assertEqual({int(r["bus"]) for r in rows if r["case"] == case}, set(range(30, 40)))

    def test_fractional_rating_and_current_semantics(self) -> None:
        with (REPORT / "tables" / "TABLE_M14_fractional_rho_semantics.csv").open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 10)
        self.assertTrue(all(r["status"] == "EVALUATED" for r in rows))
        base = {36: 700.0, 38: 1000.0}
        for row in rows:
            rho = float(row["rho"])
            self.assertTrue(math.isclose(float(row["SG_Sn_MVA"]), base[int(row["bus"])] * (1-rho), abs_tol=1e-8))
            self.assertTrue(math.isclose(float(row["GFL_port_scale"]), rho, abs_tol=1e-12))
            self.assertTrue(math.isclose(float(row["SG_P_MW"])+float(row["GFL_P_MW"]),
                                          float(row["network_P_MW"]), abs_tol=1e-8))

    def test_direct_closure_and_device_power_balance(self) -> None:
        with (REPORT / "tables" / "TABLE_M21_corrected_closure_identity.csv").open(newline="", encoding="utf-8") as fh:
            identity = list(csv.DictReader(fh))
        self.assertEqual(len(identity), 3)
        self.assertTrue(all(r["matrix_identity_pass"] == "True" and
                            r["pole_identity_pass"] == "True" for r in identity))
        with (REPORT / "tables" / "TABLE_M04_direct_device_PQ_balance.csv").open(newline="", encoding="utf-8") as fh:
            sharing = list(csv.DictReader(fh))
        self.assertEqual(len(sharing), 30)
        self.assertLess(max(abs(float(r["p_balance_MW"])) for r in sharing), 1e-8)
        self.assertLess(max(abs(float(r["q_balance_Mvar"])) for r in sharing), 1e-8)


if __name__ == "__main__":
    unittest.main()
