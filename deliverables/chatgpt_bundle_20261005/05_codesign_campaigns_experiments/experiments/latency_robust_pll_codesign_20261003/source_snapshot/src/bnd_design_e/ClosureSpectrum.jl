module ClosureSpectrum

using LinearAlgebra
using ..CollectiveModel

export port_closure, pole_sensitivities, physical_spectrum

function transfer_and_slope(J,s)
    R=(s*I(size(J.A,1))-J.A)\I(size(J.A,1))
    RB=R*J.B
    return (Y=J.D+J.C*RB,Ys=-J.C*(R*RB),R=R,RB=RB)
end

"The static 78-coordinate network is the base; every dynamic port is retained."
function port_closure(net,s,rho,Kp,Ki;derivatives=false)
    buses=sort(collect(keys(net.sg))); n=length(buses)
    length(rho)==length(Kp)==length(Ki)==n || throw(DimensionMismatch("one triple per generator"))
    P=zeros(Float64,size(net.y_static,1),2n)
    Y=zeros(ComplexF64,2n,2n); Ys=zeros(ComplexF64,2n,2n)
    d_rho=Vector{Matrix{ComplexF64}}(undef,n)
    d_kp=Vector{Matrix{ComplexF64}}(undef,n)
    d_ki=Vector{Matrix{ComplexF64}}(undef,n)
    for (j,b) in enumerate(buses)
        yi=(2b-1):(2b);pi=(2j-1):(2j);P[yi,pi].=Matrix{Float64}(I,2,2)
        sg=transfer_and_slope(net.sg[b].J,s)
        op=net.gfl[b].op
        J=CollectiveModel.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=Kp[j],ki=Ki[j])
        gf=transfer_and_slope(J,s)
        Y[pi,pi].=(1-rho[j]).*sg.Y+rho[j].*gf.Y
        Ys[pi,pi].=(1-rho[j]).*sg.Ys+rho[j].*gf.Ys
        if derivatives
            gd=CollectiveModel.AnalyticGFLPLL.gain_derivatives(op.x,op.u,op.parameters)
            dkp=J.C*gf.R*(gd.A_kp*gf.RB+gd.B_kp)
            dki=J.C*gf.R*(gd.A_ki*gf.RB+gd.B_ki)
            d_rho[j]=zeros(ComplexF64,2n,2n);d_rho[j][pi,pi].=gf.Y-sg.Y
            d_kp[j]=zeros(ComplexF64,2n,2n);d_kp[j][pi,pi].=rho[j].*dkp
            d_ki[j]=zeros(ComplexF64,2n,2n);d_ki[j][pi,pi].=rho[j].*dki
        end
    end
    G=transpose(P)*(net.y_static\P)
    C=Matrix{ComplexF64}(I,2n,2n)+Y*G
    return (C=C,Cs=Ys*G,G=G,Y=Y,P=P,
        d_rho=derivatives ? [X*G for X in d_rho] : Matrix{ComplexF64}[],
        d_kp=derivatives ? [X*G for X in d_kp] : Matrix{ComplexF64}[],
        d_ki=derivatives ? [X*G for X in d_ki] : Matrix{ComplexF64}[],
        network_condition=cond(net.y_static))
end

"Exact closure pole derivatives, valid for simple zeros away from device poles."
function pole_sensitivities(net,s,rho,Kp,Ki)
    x=port_closure(net,s,rho,Kp,Ki;derivatives=true)
    F=svd(x.C);u=F.U[:,end];v=F.V[:,end]
    den=dot(u,x.Cs*v)
    abs(den)>1e-12 || error("singular or near-defective closure pole")
    drho=[-dot(u,d*v)/den for d in x.d_rho]
    dkp=[-dot(u,d*v)/den for d in x.d_kp]
    dki=[-dot(u,d*v)/den for d in x.d_ki]
    return (ds_drho=drho,ds_dKp=dkp,ds_dKi=dki,
        closure_sigma_min=F.S[end],closure_slope=den,
        closure_residual=norm(x.C*v)/max(norm(x.C),eps()),
        right=v,left=u)
end

"Remove only a certified simple uniform-angle gauge pole."
function physical_spectrum(A;gauge_tolerance=1e-8)
    lambda=eigvals(A)
    idx=argmin(abs.(lambda))
    gauge=abs(lambda[idx])<=gauge_tolerance
    poles=gauge ? [lambda[i] for i in eachindex(lambda) if i!=idx] : collect(lambda)
    return (lambda=lambda,physical=poles,gauge_detected=gauge,
        gauge_pole=gauge ? lambda[idx] : nothing,
        spectral_abscissa=maximum(real.(poles)),
        critical=poles[argmax(real.(poles))])
end

end
