module IntermodalPathways

using LinearAlgebra

export graph_operator, psi_transform, effective_dissipative_operator,
       commutator_metric, candidate_chi_G, schur_gamma, pathway_matrix,
       block_pathway_matrix, pairwise_susceptibility, modal_diagonal, safeguarded_newton,
       pole_shift

graph_operator(Lambda::AbstractVector, Psihat::AbstractMatrix, s::Number) =
    s^2 * Matrix{ComplexF64}(I, length(Lambda), length(Lambda)) +
    Diagonal(Lambda) + Psihat

psi_transform(Phi::AbstractMatrix, A::AbstractMatrix) = Phi' * A * Phi

function effective_dissipative_operator(Psihat::AbstractMatrix, omega::Real)
    omega > 0 || throw(ArgumentError("omega must be positive"))
    return (Psihat - Psihat') / (2im * omega)
end

"""Exact block-Schur pathway aggregation for degenerate graph eigenspaces."""
function block_pathway_matrix(Psihat::AbstractMatrix,Lambda::AbstractVector,
                              s::Number,k::Integer,clusters)
    T=graph_operator(Lambda,Psihat,s)
    K=first(filter(cl->k in cl,clusters))
    R=filter(i->!(i in K),collect(eachindex(Lambda)))
    if isempty(R)
        return (target_cluster=K,complement=R,gamma=zeros(ComplexF64,length(K),length(K)),
            contributions=NamedTuple[],reconstruction_error=0.0,condition=1.0)
    end
    Trr=T[R,R]
    F=lu(Trr)
    B=T[K,R]
    C=T[R,K]
    X=F\Matrix{ComplexF64}(I,length(R),length(R))
    gamma=-(B*(F\C))
    rows=NamedTuple[]
    contribution_sum=zeros(ComplexF64,length(K),length(K))
    for (block_id,cluster) in enumerate(clusters)
        isempty(intersect(cluster,R)) && continue
        loc_modes=intersect(cluster,R)
        il=[findfirst(==(i),R) for i in loc_modes]
        for (other_id,other) in enumerate(clusters)
            isempty(intersect(other,R)) && continue
            loc_modes2=intersect(other,R)
            im=[findfirst(==(i),R) for i in loc_modes2]
            G=-(T[K,loc_modes]*X[il,im]*T[loc_modes2,K])
            contribution_sum .+= G
            push!(rows,(row_cluster=block_id,column_cluster=other_id,
                row_modes=join(loc_modes,","),column_modes=join(loc_modes2,","),
                block_trace=tr(G),block_frobenius=norm(G),block_dimension=length(K),
                diagonal_block=block_id==other_id))
        end
    end
    err=norm(gamma-contribution_sum)/max(norm(gamma),eps(Float64))
    return (target_cluster=K,complement=R,gamma=gamma,contributions=rows,
        reconstruction_error=err,condition=cond(Trr),sigma_min=minimum(svdvals(Trr)))
end

function commutator_metric(Lambda::AbstractVector, Psihat::AbstractMatrix)
    C = Diagonal(Lambda) * Psihat - Psihat * Diagonal(Lambda)
    return norm(C) / (norm(Lambda) * norm(Psihat) + eps(Float64)), C
end

function candidate_chi_G(Lambda::AbstractVector, Psihat::AbstractMatrix, s::Number)
    n = length(Lambda)
    TD = s^2 * Matrix{ComplexF64}(I, n, n) + Diagonal(Lambda) + Diagonal(diag(Psihat))
    C = Psihat - Diagonal(diag(Psihat))
    return (value=opnorm(TD \ C), sigma_min_TD=minimum(svdvals(TD)),
            condition_TD=cond(TD))
end

function schur_gamma(T::AbstractMatrix, k::Integer)
    n = size(T, 1)
    r = filter(!=(k), 1:n)
    isempty(r) && return (Tkk=T[k, k], gamma=0.0 + 0im, effective=T[k, k],
                          condition=1.0, sigma_min=Inf)
    Trr = T[r, r]
    Trk = T[r, k:k]
    Tkr = T[k:k, r]
    sigma = svdvals(Trr)
    response = Trr \ Trk
    gamma = -(Tkr * response)[1]
    eff = T[k, k] + gamma
    scalar = T[k, k] - (Tkr * response)[1]
    deterr = NaN
    if isfinite(cond(Trr)) && cond(Trr) < 1e12 && isfinite(abs(eff)) && abs(eff) > eps(Float64)
        dl, sl = logabsdet(T)
        dr, sr = logabsdet(Trr)
        rhs = dr + log(abs(eff))
        magerr = abs(dl-rhs)/max(1.0,abs(dl),abs(rhs))
        phaseerr = abs(sl-sr*(eff/abs(eff)))
        deterr = max(magerr,phaseerr)
    end
    return (Tkk=T[k, k], gamma=gamma, effective=eff,
            residual=abs(eff-scalar)/max(abs(eff),abs(scalar),eps(Float64)),
            determinant_factorization_error=deterr,
            condition=maximum(sigma)/max(minimum(sigma),eps(Float64)),
            sigma_min=minimum(sigma), complement=r)
end

function pathway_matrix(Psihat::AbstractMatrix, Lambda::AbstractVector,
                        s::Number, k::Integer)
    T = graph_operator(Lambda, Psihat, s)
    r = filter(!=(k), 1:length(Lambda))
    a = Psihat[k, r]
    b = Psihat[r, k]
    F = lu(T[r, r])
    Rb = F \ b
    # R is used only through solves: each column of the exact rank-1 pathway
    # matrix is formed from the solved response without materializing inv(Trr).
    X = F \ Matrix{ComplexF64}(I, length(r), length(r))
    G = -(reshape(a, :, 1) .* X) .* reshape(b, 1, :)
    gamma = -sum(a .* Rb)
    return (G=G, gamma=gamma, sumG=sum(G), reconstruction_error=
        abs(gamma-sum(G))/max(abs(gamma),eps(Float64)), modes=r,
        Trr=T[r,r], sigma_min=minimum(svdvals(T[r,r])), condition=cond(T[r,r]))
end

function pairwise_susceptibility(Psihat::AbstractMatrix, Lambda::AbstractVector,
                                 s::Number, k::Integer)
    rows = NamedTuple[]
    for l in eachindex(Lambda)
        l == k && continue
        tl0 = s^2 + Lambda[l] + Psihat[l, l]
        product = Psihat[k, l] * Psihat[l, k]
        push!(rows, (mode=l, coupling_product=abs(product),
                     dynamic_detuning=abs(tl0),
                     susceptibility=abs(product)/max(abs(tl0),eps(Float64)),
                     complex_susceptibility=-product/tl0, t_l0=tl0))
    end
    return sort!(rows; by=x -> -x.susceptibility)
end

modal_diagonal(Lambda::AbstractVector, Psihat::AbstractMatrix, s::Number, k::Integer) =
    s^2 + Lambda[k] + Psihat[k, k]

function safeguarded_newton(f, df, s0::Number; tol=1e-10, maxiter=80,
                           max_step=2.0, admissible=s -> true)
    s = ComplexF64(s0)
    residual = try abs(f(s)) catch; Inf end
    for iteration in 1:maxiter
        residual <= tol * max(1.0, abs(s)^2) &&
            return (root=s, converged=true, iterations=iteration-1, residual=residual)
        d = df(s)
        (!isfinite(abs(d)) || abs(d) <= 100eps(Float64)) && break
        step = f(s) / d
        abs(step) > max_step && (step *= max_step/abs(step))
        accepted = false
        alpha = 1.0
        for _ in 1:24
            trial = s - alpha*step
            if admissible(trial)
                rt = try abs(f(trial)) catch; Inf end
                if isfinite(rt) && (rt < residual || rt <= tol*max(1.0,abs(trial)^2))
                    s, residual, accepted = trial, rt, true
                    break
                end
            end
            alpha *= 0.5
        end
        accepted || break
    end
    return (root=s, converged=residual <= tol*max(1.0,abs(s)^2),
            iterations=maxiter, residual=residual)
end

function pole_shift(gamma::Number, derivative::Number, s0::Number)
    abs(derivative) > eps(Float64) || throw(ArgumentError("diagonal root derivative is singular"))
    delta = -gamma / derivative
    return (delta=delta, predicted=s0+delta)
end

end
