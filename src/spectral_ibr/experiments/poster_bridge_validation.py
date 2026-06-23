from spectral_ibr.experiments.base_gate import BaseGate, GateResult, GateRunConfig


class PosterBridgeValidationGate(BaseGate):
    name = "poster-bridge-validation"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(
            self.name,
            "not_implemented",
            message="Consumes Stage 1 processed cases and renders the pole-plane figure.",
        )
