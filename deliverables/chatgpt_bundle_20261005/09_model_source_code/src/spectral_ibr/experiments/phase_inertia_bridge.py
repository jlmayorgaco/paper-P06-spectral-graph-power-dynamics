from spectral_ibr.experiments.base_gate import GateResult, GateRunConfig


class PhaseInertiaBridgeGate:
    name = "phase_inertia_bridge"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(self.name, "not_implemented")
