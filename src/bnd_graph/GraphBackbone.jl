module GraphBackbone

using LinearAlgebra

export mass_audit, symmetric_sqrt, gauge_projector, primary_backbone,
       euclidean_backbone, generalized_basis, degenerate_clusters,
       principal_angles, mode_overlap

function mass_audit(M::AbstractMatrix; symmetry_tol=1e-12)
    residual = norm(M - M', Inf) / max(norm(M, Inf), eps(Float64))
    H = Hermitian((M + M') / 2)
    vals = eigvals(H)
    sv = svdvals(M)
    return (symmetry_error=residual, eigenvalues=vals, lambda_min=minimum(vals),
            sigma_min=minimum(sv), sigma_max=maximum(sv),
            condition=maximum(sv) / max(minimum(sv), eps(Float64)),
            spd=(residual <= symmetry_tol && minimum(vals) > 0))
end

function symmetric_sqrt(A::AbstractMatrix; inverse=false)
    E = eigen(Hermitian((A + A') / 2))
    minimum(E.values) > 0 || throw(ArgumentError("matrix must be positive definite"))
    p = inverse ? 1 ./ sqrt.(E.values) : sqrt.(E.values)
    return E.vectors * Diagonal(p) * E.vectors'
end

function gauge_projector(Mhalf::AbstractMatrix, g::AbstractVector)
    gt = Mhalf * g
    gt ./= norm(gt)
    return Matrix{Float64}(I, length(g), length(g)) - gt * gt', gt
end

"""Preregistered M-normalized Hermitian, gauge-preserving primary backbone."""
function primary_backbone(M::AbstractMatrix, L0::AbstractMatrix, g::AbstractVector)
    Mh = symmetric_sqrt(M)
    Mih = symmetric_sqrt(M; inverse=true)
    Ltilde0 = Mih * L0 * Mih
    H0 = (Ltilde0 + Ltilde0') / 2
    Pg, gt = gauge_projector(Mh, g)
    LtG = Pg * H0 * Pg
    LG = Mh * LtG * Mh
    neg_tol=1e-10*max(1.0,maximum(abs,eigvals(Hermitian(LtG));init=0.0))
    return (LG=Matrix{Float64}(real.(LG)), Ltilde=Matrix{Float64}(real.(LtG)),
            Ltilde0=Matrix{Float64}(real.(Ltilde0)), H0=Matrix{Float64}(real.(H0)),
            Mhalf=Mh, Minvhalf=Mih, projector=Pg, gauge_mass=gt,
            gauge_residual=norm(LG * g) / max(norm(LG) * norm(g), eps(Float64)),
            negative_modes=count(x -> x < -neg_tol, eigvals(Hermitian(LtG))))
end

"""Secondary Euclidean gauge-preserving Hermitian part of L0."""
function euclidean_backbone(L0::AbstractMatrix, g::AbstractVector)
    P = Matrix{Float64}(I, length(g), length(g)) - g * g' / real(dot(g, g))
    H = (L0 + L0') / 2
    LG = P * H * P
    ev=eigvals(Hermitian(LG))
    neg_tol=1e-10*max(1.0,maximum(abs,ev;init=0.0))
    return (LG=Matrix{Float64}(real.(LG)), projector=P,
            gauge_residual=norm(LG * g) / max(norm(LG) * norm(g), eps(Float64)),
            negative_modes=count(x -> x < -neg_tol, ev))
end

function degenerate_clusters(values::AbstractVector; tol=1e-8)
    groups = Vector{Vector{Int}}()
    isempty(values) && return groups
    order = sortperm(real.(values))
    current = [order[1]]
    for idx in order[2:end]
        prev = current[end]
        scale = max(1.0, abs(values[prev]), abs(values[idx]))
        if abs(values[idx] - values[prev]) <= tol * scale
            push!(current, idx)
        else
            push!(groups, current)
            current = [idx]
        end
    end
    push!(groups, current)
    return groups
end

"""Generalized basis valid for positive, zero, and negative backbone modes."""
function generalized_basis(M::AbstractMatrix, LG::AbstractMatrix; tol_deg=1e-8)
    Mh = symmetric_sqrt(M)
    Mih = symmetric_sqrt(M; inverse=true)
    Lt = Mih * LG * Mih
    Lt = (Lt + Lt') / 2
    E = eigen(Hermitian(Lt))
    p = sortperm(real.(E.values))
    Lambda = Float64.(real.(E.values[p]))
    U = Matrix{Float64}(real.(E.vectors[:, p]))
    Phi = Mih * U
    n = length(Lambda)
    eM = norm(Phi' * M * Phi - I(n), Inf)
    eL = norm(Phi' * LG * Phi - Diagonal(Lambda), Inf)
    zero_tol = 1e-10 * max(1.0, maximum(abs, Lambda; init=0.0))
    zeros = findall(abs.(Lambda) .<= zero_tol)
    neg_tol=1e-10*max(1.0,maximum(abs,Lambda;init=0.0))
    return (Phi=Phi, Lambda=Lambda, mass_modes=U, Mhalf=Mh, Minvhalf=Mih,
            M_orthogonality_error=eM, L_diagonalization_error=eL,
            zero_modes=zeros, negative_modes=findall(x -> x < -neg_tol, Lambda),
            clusters=degenerate_clusters(Lambda; tol=tol_deg))
end

function principal_angles(A::AbstractMatrix, B::AbstractMatrix)
    QA = Matrix(qr(A).Q)[:, 1:size(A, 2)]
    QB = Matrix(qr(B).Q)[:, 1:size(B, 2)]
    vals = clamp.(svdvals(QA' * QB), 0.0, 1.0)
    return acos.(vals)
end

function mode_overlap(PhiA::AbstractMatrix, MA::AbstractMatrix,
                      PhiB::AbstractMatrix, MB::AbstractMatrix)
    Mref = (MA + MB) / 2
    C = PhiA' * Mref * PhiB
    na = sqrt.(max.(real.(diag(PhiA' * Mref * PhiA)), eps(Float64)))
    nb = sqrt.(max.(real.(diag(PhiB' * Mref * PhiB)), eps(Float64)))
    return abs.(C) ./ (na * nb')
end

end
