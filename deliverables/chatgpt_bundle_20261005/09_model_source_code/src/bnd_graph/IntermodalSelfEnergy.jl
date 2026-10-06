module IntermodalSelfEnergy

using LinearAlgebra
using ..GraphModalOperator: split_sigma

export schur_self_energy, determinant_factorization_error,
       pairwise_second_order, pairwise_contributions

"""Compute the exact scalar Schur complement for graph mode k using solves."""
function schur_self_energy(T::AbstractMatrix, k::Integer)
    n = size(T, 1)
    size(T, 2) == n || throw(ArgumentError("T must be square"))
    1 <= k <= n || throw(BoundsError(T, k))
    r = filter(!=(k), collect(1:n))
    Tkk = T[k, k]
    isempty(r) && return (Tkk=Tkk, gamma=0.0 + 0.0im,
                          effective=Tkk, direct_effective=Tkk,
                          residual=0.0, complement_condition=1.0)
    Tkr = T[k:k, r]
    Trk = T[r, k:k]
    Trr = T[r, r]
    response = -(Trr \ Trk)
    gamma = (Tkr * response)[1]
    effective = Tkk + gamma
    direct_effective = Tkk + (Tkr * response)[1]
    residual = abs(effective - direct_effective) /
               max(1.0, abs(effective), abs(direct_effective))
    return (Tkk=Tkk, gamma=gamma, effective=effective,
            direct_effective=direct_effective, residual=residual,
            complement_condition=cond(Trr))
end

function determinant_factorization_error(T::AbstractMatrix, k::Integer)
    n = size(T, 1)
    r = filter(!=(k), collect(1:n))
    isempty(r) && return 0.0
    Trr = T[r, r]
    schur = schur_self_energy(T, k)
    dl, sl = logabsdet(T)
    dr, sr = logabsdet(Trr)
    rhs_log = dr + log(abs(schur.effective))
    magnitude_error = abs(dl - rhs_log) / max(1.0, abs(dl), abs(rhs_log))
    phase_rhs = sr * (schur.effective / abs(schur.effective))
    phase_error = abs(sl - phase_rhs)
    return max(magnitude_error, phase_error)
end

"""Second-order pairwise terms from the uncoupled diagonal complement."""
function pairwise_contributions(S::AbstractMatrix, Lambda::AbstractVector,
                                s::Number, k::Integer)
    parts = split_sigma(S)
    rows = NamedTuple[]
    for l in eachindex(Lambda)
        l == k && continue
        tl = s^2 + s * S[l, l] + Lambda[l]
        product = parts.offdiagonal[k, l] * parts.offdiagonal[l, k]
        gamma = abs(tl) <= eps(Float64) ? complex(NaN, NaN) :
                -s^2 * product / tl
        push!(rows, (mode_l=l, gamma=gamma, complementary_abs=abs(tl),
                     coupling_product_abs=abs(product)))
    end
    sort!(rows; by=x -> isfinite(abs(x.gamma)) ? -abs(x.gamma) : Inf)
    return rows
end

function pairwise_second_order(S::AbstractMatrix, Lambda::AbstractVector,
                               s::Number, k::Integer)
    rows = pairwise_contributions(S, Lambda, s, k)
    return sum((row.gamma for row in rows); init=0.0 + 0.0im)
end

end
