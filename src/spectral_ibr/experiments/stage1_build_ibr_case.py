from pathlib import Path

import yaml

from spectral_ibr.application.dto.ibr_case_build import IbrBuildRequest
from spectral_ibr.application.use_cases.build_ibr_case import BuildIbrCase
from spectral_ibr.experiments.base_gate import BaseGate, GateResult, GateRunConfig
from spectral_ibr.infrastructure.andes.sg_to_ibr_replacer import (
    SgToIbrReplacer,
    protocol_from_mapping,
)


class Stage1BuildIbrCaseGate(BaseGate):
    name = "stage1-build-ibr-case"

    def __init__(self, project_root: Path = Path.cwd()) -> None:
        self.project_root = project_root

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        config = config or GateRunConfig(name=self.name, case_name="mix60")
        protocol_path = config.config_path or Path("configs/cases/conversion_protocol.yaml")
        protocol_path = self._resolve(protocol_path)

        mapping = yaml.safe_load(protocol_path.read_text(encoding="utf-8"))
        protocol = protocol_from_mapping(mapping, profile_name=config.profile_name)
        case_name = config.case_name or protocol.profile_name
        output_dir = config.output_root / "stage1_build_ibr_case" / protocol.profile_name
        request = IbrBuildRequest(
            case_name=case_name,
            profile_name=protocol.profile_name,
            protocol_path=protocol_path,
            output_dir=output_dir,
            seed=config.seed,
        )
        result = BuildIbrCase(SgToIbrReplacer(self.project_root)).execute(
            request,
            protocol,
        )

        return GateResult(
            name=self.name,
            verdict=result.status,
            artifacts={key: str(path) for key, path in result.artifacts.items()},
            metrics=result.metrics,
            message=result.message,
        )

    def _resolve(self, path: Path) -> Path:
        if path.is_absolute():
            return path
        return self.project_root / path
