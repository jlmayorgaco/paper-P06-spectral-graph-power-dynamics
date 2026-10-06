module SelfEnergyCorrection

using LinearAlgebra
export gate_status, minimum_metric_correction, require_graph_seed

gate_status() = (experiment="F5",status="DIAGNOSTIC_ONLY",
    completed="Exact modal Schur self-energy Gamma_k at the frozen F2 pole",
    blocked="Gain correction and before/after performance recovery",
    reason="F4 supplies no frozen graph-controller seed theta_0 after the F2 NON_MODAL stop.")

"""Minimum-R-metric linearized residual correction from the preregistered formula."""
function minimum_metric_correction(J::AbstractMatrix,r::AbstractVector,Rinv::AbstractMatrix;rtol=1e-10)
    size(J,1)==length(r) || throw(DimensionMismatch("residual and Jacobian rows differ"))
    size(Rinv)==(size(J,2),size(J,2)) || throw(DimensionMismatch("metric size mismatch"))
    return -Rinv*J' * pinv(J*Rinv*J';rtol=rtol) * r
end

require_graph_seed() = error("F5 correction is not run: no frozen F4 graph-controller coefficient vector theta_0 exists.")

end
