from spectral_ibr.experiments.base_gate import GateResult, GateRunConfig


class PhaseBraessNepGate:
    name = "phase_braess"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(self.name, "not_implemented")
