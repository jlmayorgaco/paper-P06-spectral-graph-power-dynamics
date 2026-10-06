from spectral_ibr.experiments.base_gate import BaseGate, GateResult, GateRunConfig


class PosterWeakElementsGate(BaseGate):
    name = "poster-weak-elements"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(
            self.name,
            "not_implemented",
            message="Consumes Stage 1 processed cases and emits typed weak-element records.",
        )
