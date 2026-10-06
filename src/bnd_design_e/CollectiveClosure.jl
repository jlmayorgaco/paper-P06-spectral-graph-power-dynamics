module CollectiveClosure

using LinearAlgebra
include(joinpath(@__DIR__,"PLLLowRank.jl"))

export closure_at, determinant_audit

"""Twenty-dimensional exact multi-bus port closure away from base-device poles."""
function closure_at(net,s,rho,Kp,Ki)
    buses=sort(collect(keys(net.sg))); n=length(buses)
    length(rho)==length(Kp)==length(Ki)==n || throw(DimensionMismatch("one triple per generator"))
    ny=size(net.y_static,1); m=2n
    P=zeros(Float64,ny,m); Delta=zeros(ComplexF64,m,m)
    T0=ComplexF64.(net.y_static)
    for (j,b) in enumerate(buses)
        yi=(2b-1):(2b); pi=(2j-1):(2j)
        P[yi,pi].=I(2)
        Js=net.sg[b].J
        Ys=Js.D+Js.C*((s*I(size(Js.A,1))-Js.A)\Js.B)
        op=net.gfl[b].op
        f=PLLLowRank.pll_factors(op)
        Yg=PLLLowRank.woodbury_port(f,s,Kp[j],Ki[j])
        T0[yi,yi].+=Ys
        Delta[pi,pi].=rho[j].*(Yg-Ys)
    end
    G=transpose(P)*(T0\P)
    C=Matrix{ComplexF64}(I,m,m)+Delta*G
    T=T0+P*Delta*transpose(P)
    return (C=C,T=T,T0=T0,P=P,Delta=Delta,G=G,
        closure_dimension=m,network_dimension=ny)
end

function determinant_audit(x)
    lt,pt=logabsdet(x.T);lb,pb=logabsdet(x.T0);lc,pc=logabsdet(x.C)
    return (logabsdet_error=abs(lt-lb-lc),phase_error=abs(pt-pb*pc),
        reconstruction_error=norm(x.T-(x.T0+x.P*x.Delta*transpose(x.P)))/max(norm(x.T),eps()),
        sigma_C=minimum(svdvals(x.C)),sigma_T=minimum(svdvals(x.T)))
end

end
