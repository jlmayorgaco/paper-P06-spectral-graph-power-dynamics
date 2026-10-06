module RhoDesign

using LinearAlgebra

export rho_determinant_coefficients, rho_candidates

function subsets_of_size(n,k)
    k==0 && return [Int[]]
    k>n && return Vector{Vector{Int}}()
    out=Vector{Vector{Int}}()
    function visit(start,acc)
        if length(acc)==k; push!(out,copy(acc)); return; end
        for j in start:(n-(k-length(acc))+1)
            push!(acc,j); visit(j+1,acc); pop!(acc)
        end
    end
    visit(1,Int[])
    return out
end

"Coefficients of det(I+rho*Q), from principal minors (exact degree <= rank(Q))."
function rho_determinant_coefficients(Q::AbstractMatrix;rtol=1e-11)
    n,m=size(Q); n==m || throw(DimensionMismatch("Q must be square"))
    s=svdvals(Q); rank=count(>(rtol*maximum(s;init=0.0)),s)
    coeff=zeros(Float64,rank+1); coeff[1]=1
    for k in 1:rank
        coeff[k+1]=sum(det(Q[idx,idx]) for idx in subsets_of_size(n,k))
    end
    return coeff
end

function rho_candidates(coeff;tol=1e-10)
    length(coeff)<=1 && return Float64[]
    # Companion roots are used only on the already-derived low-degree polynomial.
    c=Float64.(coeff); scale=maximum(abs,c;init=0.0)
    while length(c)>1 && abs(last(c))<=tol*scale; pop!(c); end
    degree=length(c)-1; degree==0 && return Float64[]
    C=zeros(Float64,degree,degree)
    degree>1 && (C[2:end,1:end-1].=I(degree-1))
    C[:,end].=-c[1:degree]./c[end]
    roots=eigvals(C)
    return sort([real(z) for z in roots if abs(imag(z))<=tol*max(1,abs(real(z)))])
end

end
