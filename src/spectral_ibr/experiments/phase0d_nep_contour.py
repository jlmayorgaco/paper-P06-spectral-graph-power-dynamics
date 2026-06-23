from spectral_ibr.experiments.base_gate import GateResult, GateRunConfig


class Phase0DGate:
    name = "phase0d"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(self.name, "not_implemented")
