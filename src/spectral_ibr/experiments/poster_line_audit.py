from spectral_ibr.experiments.base_gate import BaseGate, GateResult, GateRunConfig


class PosterLineAuditGate(BaseGate):
    name = "poster-line-audit"

    def run(self, config: GateRunConfig | None = None) -> GateResult:
        return GateResult(
            self.name,
            "not_implemented",
            message="Consumes Stage 1 processed cases and emits resonant line-audit records.",
        )
