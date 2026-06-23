from spectral_ibr.experiments.base_gate import Gate
from spectral_ibr.experiments.poster_bridge_validation import PosterBridgeValidationGate
from spectral_ibr.experiments.poster_line_audit import PosterLineAuditGate
from spectral_ibr.experiments.poster_weak_elements import PosterWeakElementsGate
from spectral_ibr.experiments.phase0d_nep_contour import Phase0DGate
from spectral_ibr.experiments.phase_braess_nep import PhaseBraessNepGate
from spectral_ibr.experiments.phase_inertia_bridge import PhaseInertiaBridgeGate
from spectral_ibr.experiments.phase_passivity_defect import PhasePassivityDefectGate
from spectral_ibr.experiments.stage1_build_ibr_case import Stage1BuildIbrCaseGate
from spectral_ibr.experiments.stage1_td_compare import Stage1TdCompareGate

_GATES: dict[str, Gate] = {
    "phase0d": Phase0DGate(),
    "phase_braess": PhaseBraessNepGate(),
    "phase_inertia_bridge": PhaseInertiaBridgeGate(),
    "phase_passivity_defect": PhasePassivityDefectGate(),
    "poster-bridge-validation": PosterBridgeValidationGate(),
    "poster-line-audit": PosterLineAuditGate(),
    "poster-weak-elements": PosterWeakElementsGate(),
    "stage1-build-ibr-case": Stage1BuildIbrCaseGate(),
    "stage1-td-compare": Stage1TdCompareGate(),
}


def list_gates() -> list[str]:
    return sorted(_GATES)


def get_gate(name: str) -> Gate:
    return _GATES[name]
