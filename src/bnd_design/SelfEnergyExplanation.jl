module SelfEnergyExplanation

using LinearAlgebra

export cancellation_factor, pathway_matrix

function cancellation_factor(Gamma, pathways)
    return sum(abs,pathways)/max(abs(Gamma),eps(Float64))
end

function pathway_matrix(T::AbstractMatrix,k::Integer)
    n=size(T,1); 1<=k<=n || throw(BoundsError(T,k))
    r=[j for j in 1:n if j!=k]
    block=T[r,r]
    W=block\T[r,k]
    R=block\Matrix{ComplexF64}(I,length(r),length(r))
    g=-T[k,r]*W
    terms=zeros(ComplexF64,n,n)
    for (a,l) in enumerate(r), (b,m) in enumerate(r)
        terms[l,m]=-T[k,l]*R[a,b]*T[m,k]
    end
    return (gamma=g,pathways=terms,cancellation=cancellation_factor(g,terms))
end

end
