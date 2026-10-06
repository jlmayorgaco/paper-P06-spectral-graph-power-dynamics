using LinearAlgebra
using NetworkDynamics

export StabilityAudit, stability_audit

struct StabilityAudit
    eigenvalues::Vector{ComplexF64}
    nontrivial_eigenvalues::Vector{ComplexF64}
    gauge_eigenvalues::Vector{ComplexF64}
    max_real::Float64
    dynamic_margin::Float64
    stable::Bool
    finite::Bool
end

"""
Linearize a qualified equilibrium and report the non-gauge spectral margin.

PowerDynamics/NetworkDynamics can return a numerical zero mode at machine
precision. Such modes are recorded, not silently discarded from the raw
spectrum, and are excluded from the engineering margin only when their
absolute value is <= `gauge_tol`.
"""
function stability_audit(s0; gauge_tol = 1e-8, margin_tol = 1e-8)
    λ = ComplexF64.(jacobian_eigenvals(s0))
    finite = all(isfinite, real.(λ)) && all(isfinite, imag.(λ))
    gauge = λ[abs.(λ) .<= gauge_tol]
    nontrivial = λ[abs.(λ) .> gauge_tol]
    isempty(nontrivial) && throw(ArgumentError("no non-gauge eigenvalues remain"))
    max_real = maximum(real.(nontrivial))
    margin = -max_real
    return StabilityAudit(
        λ,
        nontrivial,
        gauge,
        max_real,
        margin,
        finite && margin > margin_tol,
        finite,
    )
end
