module DelayCharacteristic

using LinearAlgebra, ForwardDiff
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE
const R=ReducedDAE
const TAU_LPF=1/(2pi*300)

"Build the exact full-support linear DDE characteristic data at one frozen design."
function linearization(ctx,rho,kp,ki)
    m=R.model(ctx,rho,kp,ki;bus=8,delta=0.,dc_convention=:physical_supply)
    s=R.N.spectrum(ctx,rho,kp,ki)
    Q=s.quotient; A=Q' * s.model.Ared * Q
    Ai=Matrix{Float64}[]; bcols=Vector{Float64}[]; bpcols=Vector{Float64}[]; bicolumns=Vector{Float64}[]; ccols=Vector{Float64}[]; labels=Int[]
    for i in eachindex(rho)
        rho[i]>0 || continue
        ix=m.gfidx[i]; busport=(2i-1):(2i)
        c=ForwardDiff.gradient(m.x0) do x
            v=R.voltage(x,m);theta=x[ix[3]]
            -sin(theta)*v[first(busport)]+cos(theta)*v[last(busport)]
        end
        bp=zeros(length(m.x0));bp[ix[4]]=1/TAU_LPF
        bi=zeros(length(m.x0));bi[ix[5]]=1
        b=kp[i].*bp+ki[i].*bi
        bq=Q'*b;cq=Q'*c
        push!(bcols,bq);push!(bpcols,Q'*bp);push!(bicolumns,Q'*bi);push!(ccols,cq);push!(Ai,bq*cq');push!(labels,29+i)
    end
    A0=A-sum(Ai;init=zeros(size(A)))
    (;m,Q,A,A0,Ai,B=hcat(bcols...),Bp=hcat(bpcols...),Bi=hcat(bicolumns...),C=hcat(ccols...),labels,
      alpha=maximum(real,eigvals(A)),rho=copy(rho),kp=copy(kp),ki=copy(ki))
end

"Exact retarded-DDE characteristic matrix Delta(s)=sI-A0-sum(Ai*exp(-s*taui))."
function delta_matrix(L,s,tau)
    length(tau)==length(L.Ai) || throw(DimensionMismatch("one fixed delay per active GFL PLL channel"))
    T=promote_type(typeof(s),ComplexF64)
    D=Matrix{T}(s*I-L.A0)
    for j in eachindex(L.Ai)
        D .-= exp(-s*tau[j]).*L.Ai[j]
    end
    D
end

function delta_s(L,s,tau)
    D=Matrix{promote_type(typeof(s),ComplexF64)}(I,size(L.A,1),size(L.A,2))
    for j in eachindex(L.Ai)
        D .+= (tau[j]*exp(-s*tau[j])).*L.Ai[j]
    end
    D
end

"Exact derivative of Delta with respect to one fixed exogenous delay."
delta_tau(L,s,tau,j) = s*exp(-s*tau[j]).*L.Ai[j]

"Normalized nonlinear-eigenvalue residual for a candidate right vector."
function residual(L,s,tau,v)
    D=delta_matrix(L,s,tau)
    norm(D*v)/(max(1,norm(D)*norm(v)))
end

end
