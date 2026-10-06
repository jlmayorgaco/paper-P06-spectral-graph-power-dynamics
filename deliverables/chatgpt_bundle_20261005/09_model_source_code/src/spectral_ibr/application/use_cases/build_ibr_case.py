from spectral_ibr.application.dto.ibr_case_build import (
    IbrBuildRequest,
    IbrBuildResult,
    IbrConversionProtocol,
)
from spectral_ibr.application.ports.ibr_case_replacer import IIbrCaseReplacer


class BuildIbrCase:
    """Stage 1 use case: turn a citable conversion protocol into artifacts."""

    def __init__(self, replacer: IIbrCaseReplacer) -> None:
        self.replacer = replacer

    def execute(
        self,
        request: IbrBuildRequest,
        protocol: IbrConversionProtocol,
    ) -> IbrBuildResult:
        return self.replacer.build(request, protocol)
