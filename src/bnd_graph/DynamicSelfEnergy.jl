module DynamicSelfEnergy

using LinearAlgebra

export StateSpaceSelfEnergy, sigma, sigma_derivative, augmented_matrix,
       diagonal_only, with_offdiagonal_scale

"""Real state-space realization Σ(s) = D₀ + C(sI-Ac)⁻¹B.

`channel_pairs[r] = (i,j)` records that controller state r realizes a
transfer path from retained input j to retained output i. It enables a
well-defined diagonal modal projection for the synthetic realizations.
"""
struct StateSpaceSelfEnergy
    D0::Matrix{Float64}
    Ac::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    channel_pairs::Vector{Tuple{Int,Int}}
    offdiag_state_indices::Vector{Int}
    offdiag_base_scale::Float64
end

function StateSpaceSelfEnergy(D0::AbstractMatrix, Ac::AbstractMatrix,
                              B::AbstractMatrix, C::AbstractMatrix;
                              channel_pairs=Tuple{Int,Int}[],
                              offdiag_state_indices=Int[],
                              offdiag_base_scale::Real=1.0)
    n = size(D0, 1)
    size(D0) == (n, n) || throw(ArgumentError("D₀ must be square"))
    size(Ac, 1) == size(Ac, 2) || throw(ArgumentError("Ac must be square"))
    size(B) == (size(Ac, 1), n) || throw(ArgumentError("B has incompatible dimensions"))
    size(C) == (n, size(Ac, 1)) || throw(ArgumentError("C has incompatible dimensions"))
    isempty(channel_pairs) || length(channel_pairs) == size(Ac, 1) ||
        throw(ArgumentError("channel_pairs must have one entry per controller state"))
    return StateSpaceSelfEnergy(Matrix{Float64}(D0), Matrix{Float64}(Ac),
                                Matrix{Float64}(B), Matrix{Float64}(C),
                                collect(channel_pairs), collect(offdiag_state_indices),
                                Float64(offdiag_base_scale))
end

function sigma(sys::StateSpaceSelfEnergy, s::Number)
    n = size(sys.D0, 1)
    size(sys.Ac, 1) == 0 && return complex.(sys.D0)
    R = s .* Matrix{ComplexF64}(I, size(sys.Ac, 1), size(sys.Ac, 1)) .- sys.Ac
    return complex.(sys.D0) + sys.C * (R \ sys.B)
end

function sigma_derivative(sys::StateSpaceSelfEnergy, s::Number)
    n = size(sys.D0, 1)
    size(sys.Ac, 1) == 0 && return zeros(ComplexF64, n, n)
    R = s .* Matrix{ComplexF64}(I, size(sys.Ac, 1), size(sys.Ac, 1)) .- sys.Ac
    X = R \ sys.B
    return -sys.C * (R \ X)
end

function augmented_matrix(M::AbstractMatrix, L::AbstractMatrix,
                          sys::StateSpaceSelfEnergy)
    n = size(M, 1)
    size(L) == (n, n) || throw(ArgumentError("L has incompatible dimensions"))
    size(sys.D0) == (n, n) || throw(ArgumentError("self-energy has incompatible dimensions"))
    nc = size(sys.Ac, 1)
    Znn = zeros(Float64, n, n)
    Znc = zeros(Float64, n, nc)
    Zcn = zeros(Float64, nc, n)
    I_n = Matrix{Float64}(I, n, n)
    return [Znn I_n Znc;
            -(M \ L) -(M \ sys.D0) -(M \ sys.C);
            Zcn sys.B sys.Ac]
end

function diagonal_only(sys::StateSpaceSelfEnergy; include_direct_diagonal::Bool=true)
    n = size(sys.D0, 1)
    nc = size(sys.Ac, 1)
    D = include_direct_diagonal ? Matrix(Diagonal(diag(sys.D0))) : zeros(n, n)
    B = copy(sys.B)
    C = copy(sys.C)
    if isempty(sys.channel_pairs) && nc > 0
        throw(ArgumentError("diagonal projection requires channel_pairs metadata"))
    end
    for r in 1:nc
        i, j = sys.channel_pairs[r]
        if i != j
            B[r, :] .= 0.0
            C[:, r] .= 0.0
        else
            for a in 1:n
                a == j || (B[r, a] = 0.0)
                a == i || (C[a, r] = 0.0)
            end
        end
    end
    return StateSpaceSelfEnergy(D, sys.Ac, B, C;
        channel_pairs=sys.channel_pairs,
        offdiag_state_indices=Int[], offdiag_base_scale=0.0)
end

"""Return a copy with all declared off-diagonal controller paths scaled by ε."""
function with_offdiagonal_scale(sys::StateSpaceSelfEnergy, epsilon::Real)
    epsilon >= 0 || throw(ArgumentError("epsilon must be nonnegative"))
    B = copy(sys.B)
    C = copy(sys.C)
    ratio = sys.offdiag_base_scale == 0 ? 0.0 : sqrt(epsilon / sys.offdiag_base_scale)
    for r in sys.offdiag_state_indices
        B[r, :] .*= ratio
        C[:, r] .*= ratio
    end
    return StateSpaceSelfEnergy(sys.D0, sys.Ac, B, C;
        channel_pairs=sys.channel_pairs,
        offdiag_state_indices=sys.offdiag_state_indices,
        offdiag_base_scale=Float64(epsilon))
end

end
