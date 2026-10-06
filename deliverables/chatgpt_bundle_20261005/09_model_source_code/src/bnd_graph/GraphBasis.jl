module GraphBasis

using LinearAlgebra

export GraphBasisResult, generalized_graph_basis, degenerate_clusters,
       block_frobenius, oscillatory_indices

struct GraphBasisResult
    Phi::Matrix{Float64}
    Lambda::Vector{Float64}
    M_orthogonality_error::Float64
    L_diagonalization_error::Float64
    zero_modes::Vector{Int}
    clusters::Vector{Vector{Int}}
end

function degenerate_clusters(values::AbstractVector{<:Real}; tol_deg::Real=1e-8)
    n = length(values)
    groups = Vector{Vector{Int}}()
    i = 1
    while i <= n
        group = [i]
        j = i + 1
        while j <= n
            scale = max(1.0, abs(values[i]), abs(values[j]))
            abs(values[j] - values[i]) <= tol_deg * scale || break
            push!(group, j)
            j += 1
        end
        push!(groups, group)
        i = j
    end
    return groups
end

"""Mass-normalized generalized eigenbasis for L φ = ν M φ.

The method uses a Cholesky factor M = R'R and never forms M⁻¹ or M⁻¹/².
"""
function generalized_graph_basis(M::AbstractMatrix, L::AbstractMatrix;
                                 zero_tol::Real=1e-10,
                                 tol_deg::Real=1e-8,
                                 symmetry_tol::Real=1e-11)
    size(M, 1) == size(M, 2) || throw(ArgumentError("M must be square"))
    size(L) == size(M) || throw(ArgumentError("L and M must have the same size"))
    norm(M - M', Inf) <= symmetry_tol * max(1.0, norm(M, Inf)) ||
        throw(ArgumentError("M must be Hermitian"))
    norm(L - L', Inf) <= symmetry_tol * max(1.0, norm(L, Inf)) ||
        throw(ArgumentError("L must be Hermitian"))

    Mf = cholesky(Hermitian((M + M') / 2))
    R = Matrix(Mf.U)
    Lmass = (R' \ Matrix(L)) / R
    Lmass = (Lmass + Lmass') / 2
    E = eigen(Hermitian(Lmass))
    order = sortperm(real.(E.values))
    values = Float64.(real.(E.values[order]))
    Q = Matrix{Float64}(real.(E.vectors[:, order]))
    scale = max(1.0, maximum(abs, values; init=0.0))
    minimum(values; init=0.0) >= -zero_tol * scale ||
        throw(ArgumentError("L must be positive semidefinite on the modeled subspace"))
    values[abs.(values) .<= zero_tol * scale] .= 0.0

    Phi = R \ Q
    Lambda = Diagonal(values)
    I_n = Matrix{Float64}(I, size(M, 1), size(M, 1))
    eM = norm(Phi' * M * Phi - I_n, Inf)
    eL = norm(Phi' * L * Phi - Lambda, Inf)
    zeros = findall(v -> abs(v) <= zero_tol * scale, values)
    return GraphBasisResult(Phi, values, eM, eL, zeros,
                            degenerate_clusters(values; tol_deg=tol_deg))
end

"""Frobenius norm of an operator block, invariant to unitary basis rotations
inside the supplied degenerate eigenspace.
"""
block_frobenius(A::AbstractMatrix, rows::AbstractVector{<:Integer},
                cols::AbstractVector{<:Integer}=rows) = norm(A[rows, cols])

oscillatory_indices(basis::GraphBasisResult) =
    setdiff(collect(eachindex(basis.Lambda)), basis.zero_modes)

end
