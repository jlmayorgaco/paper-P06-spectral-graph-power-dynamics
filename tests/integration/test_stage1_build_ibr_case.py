from pathlib import Path

from openpyxl import Workbook
import yaml

from spectral_ibr.application.dto.ibr_case_build import IbrBuildRequest
from spectral_ibr.application.use_cases.build_ibr_case import BuildIbrCase
from spectral_ibr.domain.errors import MissingBaseCaseError
from spectral_ibr.infrastructure.andes.sg_to_ibr_replacer import (
    SgToIbrReplacer,
    protocol_from_mapping,
)


def test_stage1_build_writes_plan_and_manifest(tmp_path: Path) -> None:
    base_case = tmp_path / "base.xlsx"
    _write_minimal_andes_case(base_case)
    protocol_path = tmp_path / "conversion_protocol.yaml"
    protocol_path.write_text(
        yaml.safe_dump(
            {
                "name": "test_protocol",
                "base_case": str(base_case),
                "active_profile": "mix60",
                "eligible_buses": [31, 32],
                "preserve_power_flow": True,
                "profiles": {
                    "mix60": {
                        "output_case": str(tmp_path / "mix60.xlsx"),
                        "keep_synchronous": [32, 39],
                        "replacements": {
                            "gfl": {
                                "gens": [2],
                                "buses": [31],
                                "model": "REGCP1",
                                "with_pll": True,
                            }
                        },
                    },
                },
                "citations": ["andes_rengen_models"],
            }
        ),
        encoding="utf-8",
    )
    protocol = protocol_from_mapping(yaml.safe_load(protocol_path.read_text()))
    request = IbrBuildRequest("mix60", "mix60", protocol_path, tmp_path / "out")

    result = BuildIbrCase(SgToIbrReplacer(tmp_path)).execute(request, protocol)

    assert result.status == "built"
    assert result.artifacts["processed_case"].exists()
    assert result.artifacts["conversion_plan"].exists()
    assert result.manifest_path.exists()


def test_stage1_build_fails_without_base_case(tmp_path: Path) -> None:
    protocol = protocol_from_mapping(
        {
            "name": "test_protocol",
            "base_case": str(tmp_path / "missing.xlsx"),
            "active_profile": "mix60",
            "eligible_buses": [31],
            "preserve_power_flow": True,
            "profiles": {
                "mix60": {
                    "output_case": str(tmp_path / "mix60.xlsx"),
                    "keep_synchronous": [],
                    "replacements": {
                        "gfl": {
                            "gens": [2],
                            "buses": [31],
                            "model": "REGCP1",
                            "with_pll": True,
                        }
                    },
                },
            },
        }
    )
    request = IbrBuildRequest("mix60", "mix60", tmp_path / "protocol.yaml", tmp_path / "out")

    try:
        BuildIbrCase(SgToIbrReplacer(tmp_path)).execute(request, protocol)
    except MissingBaseCaseError as exc:
        assert "No fallback is used" in str(exc)
    else:
        raise AssertionError("Expected MissingBaseCaseError")


def _write_minimal_andes_case(path: Path) -> None:
    wb = Workbook()
    default = wb.active
    wb.remove(default)
    sheets = {
        "PV": [
            ["uid", "idx", "u", "name", "Sn", "Vn", "bus", "busr", "p0", "q0", "pmax", "pmin", "qmax", "qmin", "v0", "vmax", "vmin", "ra", "xs"],
            [0, 2, 1, 2, 100, 13.8, 31, None, 1, 0, 2, 0, 1, -1, 1.0, 1.4, 0.6, 0, 0.2],
        ],
        "Slack": [
            ["uid", "idx", "u", "name", "Sn", "Vn", "bus", "busr", "p0", "q0", "pmax", "pmin", "qmax", "qmin", "v0", "vmax", "vmin", "ra", "xs", "a0"],
            [0, 10, 1, 10, 100, 13.8, 39, None, 1, 0, 2, 0, 1, -1, 1.0, 1.4, 0.6, 0, 0.2, 0],
        ],
        "GENROU": [
            ["uid", "idx", "u", "name", "bus", "gen", "coi", "Sn", "Vn", "fn", "D", "M"],
            [0, "GENROU_2", 1, "GENROU_2", 31, 2, None, 100, 13.8, 60, 0, 5],
        ],
        "TGOV1N": [
            ["uid", "idx", "u", "name", "syn"],
            [0, "TGOV1_2", 1, "TGOV1_2", "GENROU_2"],
        ],
        "IEEEX1": [
            ["uid", "idx", "u", "name", "syn"],
            [0, "IEEEX1_2", 1, "IEEEX1_2", "GENROU_2"],
        ],
        "IEEEST": [
            ["uid", "idx", "u", "name", "avr", "busf"],
            [0, "IEEEST_2", 1, "IEEEST_2", "IEEEX1_2", "BusFreq_2"],
        ],
        "BusFreq": [
            ["uid", "idx", "u", "name", "bus", "Tf", "Tw", "fn"],
            [0, "BusFreq_2", 1, "BusFreq_2", 31, 0.02, 0.02, 60],
        ],
    }
    for name, rows in sheets.items():
        ws = wb.create_sheet(name)
        for row in rows:
            ws.append(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)
