from spectral_ibr.experiments.base_gate import GateResult, GateRunConfig


class PhasePassivityDefectGate:
    name = "phase_passivity_defect"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(self.name, "not_implemented")
