module GraphModalOperator

using LinearAlgebra

export graph_sigma, graph_sigma_derivative, graph_operator, split_sigma,
       hermitian_dissipative, antihermitian_reactive, offdiagonal_ratio,
       commutator_diagnostic, harmonic_power_identity

graph_sigma(basis, sigma_fun, s::Number) =
    basis.Phi' * sigma_fun(s) * basis.Phi

graph_sigma_derivative(basis, derivative_fun, s::Number) =
    basis.Phi' * derivative_fun(s) * basis.Phi

function graph_operator(basis, sigma_fun, s::Number)
    n = length(basis.Lambda)
    S = graph_sigma(basis, sigma_fun, s)
    return s^2 .* Matrix{ComplexF64}(I, n, n) + s .* S + Diagonal(basis.Lambda)
end

function split_sigma(S::AbstractMatrix)
    D = Matrix(Diagonal(diag(S)))
    return (diagonal=D, offdiagonal=Matrix(S - D))
end

hermitian_dissipative(S::AbstractMatrix) = (S + S') / 2
antihermitian_reactive(S::AbstractMatrix) = (S - S') / (2im)

function offdiagonal_ratio(S::AbstractMatrix; eps_guard::Real=eps(Float64))
    parts = split_sigma(S)
    return norm(parts.offdiagonal) /
           (norm(S) + eps_guard)
end

function commutator_diagnostic(Lambda::AbstractVector, S::AbstractMatrix;
                               eps_guard::Real=eps(Float64))
    L = Diagonal(Lambda)
    C = L * S - S * L
    denom = norm(L) * norm(S) + eps_guard
    return (commutator=C, norm=norm(C), chi=norm(C) / denom)
end

"""Compare direct harmonic average power with ω² zᴴHerm(Σ̂)z / 2.

The convention takes the self-energy force on the left side of the retained
equation; positive computed power is therefore power removed from the
mechanical coordinates.
"""
function harmonic_power_identity(Phi::AbstractMatrix, sigma_fun, omega::Real,
                                 z::AbstractVector)
    v = im * omega .* (Phi * z)
    f = sigma_fun(im * omega) * v
    p_direct = real(dot(v, f)) / 2
    Shat = Phi' * sigma_fun(im * omega) * Phi
    DG = hermitian_dissipative(Shat)
    p_graph = omega^2 * real(dot(z, DG * z)) / 2
    relerr = abs(p_direct - p_graph) / max(1.0, abs(p_direct), abs(p_graph))
    return (direct=p_direct, graph=p_graph, relative_error=relerr)
end

end
