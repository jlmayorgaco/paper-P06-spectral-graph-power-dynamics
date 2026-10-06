"""Evidence gates for the frozen Experiment N result."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def rows(name: str) -> list[dict]:
    with (OUT / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


class ExperimentNGates(unittest.TestCase):
    def test_model_identity_and_trim(self) -> None:
        self.assertEqual(len(rows("TABLE_N04_trim_validation.csv")), 60)
        self.assertTrue(all(r["status"] == "PASS" for r in
                            rows("TABLE_N04_trim_validation.csv")))
        mats = rows("TABLE_N05_parametric_identity.csv")
        ports = rows("TABLE_N06_withheld_port_identity.csv")
        poles = rows("TABLE_N07_withheld_pole_identity.csv")
        self.assertEqual(len(mats), 69)
        self.assertTrue(all(r["matrix_gate"] == "PASS" for r in mats))
        self.assertTrue(all(r["port_gate"] == "PASS" for r in ports))
        self.assertTrue(all(r["pole_gate"] == "PASS" for r in poles))
        self.assertLess(max(float(r["max_matched_pole_error"]) for r in poles), 1e-6)
        self.assertTrue(all(r["classification"] != "UNRESOLVED" for r in
                            rows("TABLE_N03_near_zero_mode_classification.csv")))
        grads = rows("TABLE_N08_derivative_validation.csv")
        self.assertEqual(len(grads), 90)
        self.assertTrue(all(r["status"] == "PASS" for r in grads))

    def test_frozen_candidate_and_local_certificate(self) -> None:
        path = OUT / "Z_N_NOMINAL_FINAL.toml"
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         (OUT / (path.name + ".sha256")).read_text().strip())
        candidate = tomllib.loads(path.read_text())
        model = json.loads((OUT / "MODEL_FREEZE.json").read_text())
        self.assertEqual(candidate["model_sha"], model["MODEL_SHA"])
        self.assertEqual(candidate["support_bus"], 38)
        self.assertLess(candidate["alpha_analytic_per_s"], -0.05)
        self.assertFalse(candidate["global_certified"])
        self.assertEqual(len(rows("TABLE_N09_support_lattice_budget_screen.csv")), 2048)
        self.assertTrue(all(r["feasible"].lower() == "false" for r in
                            rows("TABLE_N09_support_lattice_budget_screen.csv")))
        self.assertTrue(any(r["KKT"].lower() == "true" and int(r["bus"]) == 38
                            for r in rows("TABLE_N10_all_KKT_candidates.csv")))

    def test_postfreeze_independent_validation(self) -> None:
        pd = rows("TABLE_N13_powerdynamics_postfreeze_validation.csv")[0]
        self.assertEqual(pd["PD_validation"], "PASS")
        self.assertLess(float(pd["PD_alpha"]), -0.05)
        self.assertLess(float(pd["complete_finite_pole_max_error"]), 1e-5)
        tds = rows("TABLE_N14_time_domain_validation.csv")
        self.assertEqual({int(r["event_bus"]) for r in tds}, {8, 16, 29})
        self.assertEqual(len(tds), 6)
        self.assertTrue(all(r["TDS_validation"] == "PASS" for r in tds))
        self.assertLess(max(float(r["frequency_scaling_error"]) for r in tds), .05)
        self.assertTrue(math.isfinite(float(pd["critical_right_eigenvector_overlap"])))


if __name__ == "__main__":
    unittest.main()
