module LowRankUpdates

using LinearAlgebra
using ..AnalyticDeviceModel: rank_audit

export exact_svd_factor, affine_reconstruction, woodbury_resolvent, woodbury_error

"Numerical SVD factor with full-spectrum rank and reconstruction diagnostics."
function exact_svd_factor(A::AbstractMatrix; rtol=nothing)
    F=svd(Matrix{Float64}(A))
    audit=rank_audit(A;rtol)
    r=audit.rank
    U=F.U[:,1:r]*Diagonal(F.S[1:r])
    V=F.V[:,1:r]
    return (U=U,V=V,rank=r,relative_residual=norm(A-U*V')/max(norm(A),eps(Float64)),
            singular_values=F.S,rank_tolerance=audit.tolerance)
end

function affine_reconstruction(A0,updates,parameters)
    size(updates,3)==length(parameters) || throw(DimensionMismatch("one parameter per update slice required"))
    A=copy(A0)
    for j in eachindex(parameters); A .+= parameters[j].*view(updates,:,:,j); end
    return A
end

function woodbury_resolvent(A0::AbstractMatrix,U::AbstractMatrix,Theta::AbstractMatrix,
                            V::AbstractMatrix,s::Number)
    n=size(A0,1); size(A0,2)==n || throw(DimensionMismatch("A0 must be square"))
    size(U,1)==n && size(V,1)==n || throw(DimensionMismatch("update factors have wrong row count"))
    size(U,2)==size(V,2)==size(Theta,1)==size(Theta,2) ||
        throw(DimensionMismatch("Woodbury update dimensions disagree"))
    F=lu(s*I-A0)
    R0=F\Matrix{ComplexF64}(I,n,n)
    R0U=F\ComplexF64.(U)
    W=I-Theta*(V'*R0U)
    return R0+R0U*(W\(Theta*(V'*R0)))
end

function woodbury_error(A0,U,Theta,V,s)
    A=A0+U*Theta*V'
    Rd=inv(s*I-A)
    Rw=woodbury_resolvent(A0,U,Theta,V,s)
    return norm(Rd-Rw)/max(norm(Rd),eps(Float64))
end

end
