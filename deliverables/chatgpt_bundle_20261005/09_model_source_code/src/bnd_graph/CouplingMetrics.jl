module CouplingMetrics

using LinearAlgebra
using ..GraphModalOperator: graph_sigma, hermitian_dissipative,
                            offdiagonal_ratio, commutator_diagnostic

export graph_frequency_metrics, candidate_chiG, frequency_grid_hz,
       spectral_abscissa

function frequency_grid_hz(config=nothing)
    cfg = config === nothing ? Dict{String,Any}() : config
    fmin = Float64(get(cfg,"minimum",0.01))
    fmax = Float64(get(cfg,"maximum",100.0))
    dmin = Float64(get(cfg,"dense_minimum",0.1))
    dmax = Float64(get(cfg,"dense_maximum",10.0))
    coarse = 10.0 .^ range(log10(fmin), log10(fmax); length=120)
    dense = 10.0 .^ range(log10(dmin), log10(dmax); length=360)
    return sort!(unique!(vcat(collect(coarse), collect(dense))))
end

function graph_frequency_metrics(Lambda::AbstractVector, S::AbstractMatrix)
    DG = hermitian_dissipative(S)
    comm = commutator_diagnostic(Lambda, S)
    return (sigma_norm2=opnorm(S, 2), sigma_normF=norm(S),
            offdiag_normF=norm(S - Diagonal(diag(S))),
            offdiag_ratio=offdiagonal_ratio(S),
            commutator_normF=comm.norm, chi_comm=comm.chi,
            DG_lambda_min=minimum(real.(eigvals(Hermitian(DG)))),
            DG_lambda_max=maximum(real.(eigvals(Hermitian(DG)))),
            DG=DG)
end

"""Axis-sampled candidate metric ||T_D(jω)⁻¹ C_G(jω)||₂.

This finite-grid value is a diagnostic, not a stability theorem.
"""
function candidate_chiG(Lambda::AbstractVector, sigma_fun, frequencies_hz::AbstractVector)
    n = length(Lambda)
    rows = NamedTuple[]
    for f in frequencies_hz
        omega = 2π * f
        s = im * omega
        S = sigma_fun(s)
        D = Matrix(Diagonal(diag(S)))
        OD = S - D
        TD = s^2 .* Matrix{ComplexF64}(I, n, n) + s .* D + Diagonal(Lambda)
        CG = s .* OD
        value = try
            opnorm(TD \ CG, 2)
        catch
            Inf
        end
        push!(rows, (frequency_hz=f, chiG=value))
    end
    idx = argmax(getfield.(rows, :chiG))
    return (samples=rows, supremum=rows[idx].chiG,
            frequency_at_supremum_hz=rows[idx].frequency_hz)
end

spectral_abscissa(A::AbstractMatrix) = maximum(real.(eigvals(A)))

end
