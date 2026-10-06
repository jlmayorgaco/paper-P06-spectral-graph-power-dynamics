module ModalGFL

using LinearAlgebra

export realify_admittance, modal_transform, offdiagonal_block_ratio,
       characteristic_pencil, exact_modal_self_energy, modal_structure

"""Interleaved rectangular-coordinate realification of a complex nodal matrix."""
function realify_admittance(Y::AbstractMatrix{<:Complex})
    n,m=size(Y); n==m || throw(DimensionMismatch("Y must be square"))
    R=zeros(ComplexF64,2n,2n)
    for i in 1:n,j in 1:n
        z=Y[i,j]
        R[2i-1,2j-1]=real(z); R[2i-1,2j]=-imag(z)
        R[2i,2j-1]=imag(z); R[2i,2j]=real(z)
    end
    return R
end

"""Node-mode transform that preserves interleaved real/imaginary channels."""
modal_transform(U::AbstractMatrix) = kron(U,Matrix{Float64}(I,2,2))

"""Relative Frobenius norm of inter-mode 2x2 blocks."""
function offdiagonal_block_ratio(A::AbstractMatrix; block_size=2)
    n,m=size(A); n==m || throw(DimensionMismatch("A must be square"))
    n%block_size==0 || throw(DimensionMismatch("matrix dimension must be block aligned"))
    nb=n÷block_size
    off=copy(A)
    for k in 1:nb
        ix=((k-1)*block_size+1):(k*block_size)
        off[ix,ix].=0
    end
    return norm(off)/max(norm(A),eps(Float64)),off
end

"""Assemble the exact homogeneous descriptor characteristic matrix P(s).

Positive terminal current is defined as device current entering the network,
so KCL is `Ynetwork*v + i_device = 0`. The determinant is a polynomial of
degree at most `n*nstates` in s after algebraic constraints are included.
"""
function characteristic_pencil(s,A,B,C,D,Ynetwork)
    n=size(Ynetwork,1)÷2
    nstates=size(A,1)
    size(Ynetwork)==(2n,2n) || throw(DimensionMismatch("network channel size mismatch"))
    nI=Matrix{Float64}(I,n,n)
    IA=kron(nI,Matrix(s*I(nstates)-A))
    IB=kron(nI,Matrix(B))
    IC=kron(nI,Matrix(C))
    ID=kron(nI,Matrix(D))
    top=[IA -IB]
    bottom=[IC (Ynetwork+ID)]
    return [top;bottom]
end

"""Exact Schur self-energy for each 2-channel graph-mode block."""
function exact_modal_self_energy(Tmodal::AbstractMatrix; block_size=2)
    n=size(Tmodal,1); size(Tmodal,2)==n || throw(DimensionMismatch("Tmodal must be square"))
    n%block_size==0 || throw(DimensionMismatch("Tmodal block size mismatch"))
    nb=n÷block_size; rows=NamedTuple[]
    for k in 1:nb
        keep=collect(((k-1)*block_size+1):(k*block_size))
        rest=setdiff(collect(1:n),keep)
        Tkk=Tmodal[keep,keep]
        Γ=isempty(rest) ? zeros(eltype(Tmodal),block_size,block_size) :
            -Tmodal[keep,rest]*(Tmodal[rest,rest]\Tmodal[rest,keep])
        push!(rows,(mode=k,Tkk=Matrix(Tkk),Gamma=Matrix(Γ),
            effective=Matrix(Tkk+Γ),schur_residual=norm(Tkk+Γ-
                (Tmodal[keep,keep]-Tmodal[keep,rest]*(Tmodal[rest,rest]\Tmodal[rest,keep])))))
    end
    return rows
end

"""Classify the graph basis from its transformed network coupling."""
function modal_structure(Lc,Yport; exact_tol=1e-10, approximate_tol=0.05)
    E=eigen(Hermitian((Lc+Lc')/2)); p=sortperm(E.values)
    U=Matrix(E.vectors[:,p]); G=real.(Yport); B=imag.(Yport)
    Gm=U'*G*U; Bm=U'*B*U
    ynorm=max(norm(Yport),eps(Float64))
    Goff=norm(Gm-Diagonal(diag(Gm)))/ynorm
    Boff=norm(Bm-Diagonal(diag(Bm)))/ynorm
    comm=norm(G*B-B*G)/max(norm(G)*norm(B),eps(Float64))
    classification=Goff<=exact_tol && Boff<=exact_tol ? "EXACT_SCALAR_MODAL" :
        max(Goff,Boff)<=approximate_tol ? "APPROX_MODAL_WITH_SMALL_COUPLING" : "NON_MODAL"
    return (U=U,lambda=Float64.(E.values[p]),Gmodal=Gm,Bmodal=Bm,
        conductance_offdiag_ratio=Goff,susceptance_offdiag_ratio=Boff,
        commutator_ratio=comm,classification=classification)
end

end
