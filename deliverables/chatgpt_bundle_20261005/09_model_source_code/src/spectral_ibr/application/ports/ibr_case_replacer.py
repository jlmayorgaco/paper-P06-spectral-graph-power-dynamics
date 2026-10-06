from typing import Protocol

from spectral_ibr.application.dto.ibr_case_build import (
    IbrBuildRequest,
    IbrBuildResult,
    IbrConversionProtocol,
)


class IIbrCaseReplacer(Protocol):
    def build(
        self,
        request: IbrBuildRequest,
        protocol: IbrConversionProtocol,
    ) -> IbrBuildResult:
        """Build or plan a deterministic SG-to-IBR case conversion."""
