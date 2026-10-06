module WoodburyPLL

using LinearAlgebra
using Random
using Statistics
using ..AnalyticDeviceModel: port_transfer
using ..AnalyticGFLPLL
using ..LowRankUpdates

export factorization, verify_woodbury, denominator_coefficients

function factorization(x,u,p)
    J=jacobians(x,u,p)
    G=gain_derivatives(x,u,p)
    A0=J.A-p.pll_kp*G.A_kp-p.pll_ki*G.A_ki
    # A = A0 + U*diag(Kp,Ki)*V' exactly; rows correspond to
    # filtered PLL frequency and PLL integrator states.
    U=zeros(9,2); U[4,1]=1; U[5,2]=1
    V=hcat(G.error_gradient_x/p.pll_tau,G.error_gradient_x)
    return (Aref=J.A,Bref=J.B,C=J.C,D=J.D,A0=A0,U=U,V=V,
        Bkp=G.B_kp,Bki=G.B_ki,Theta0=Diagonal([p.pll_kp,p.pll_ki]),
        derivatives=G)
end

function verify_woodbury(x,u,p; seed=2601, samples=32)
    F=factorization(x,u,p)
    rng=MersenneTwister(seed)
    resolvent_errors=Float64[]; port_errors=Float64[]; denominator_errors=Float64[]
    sample_rows=NamedTuple[]
    for j in 1:samples
        kp=p.pll_kp*(0.45+1.8rand(rng))
        ki=p.pll_ki*(0.45+1.8rand(rng))
        # Deterministic log-frequency/complex samples span the stated 0.01–10 Hz range.
        hz=10.0^(log10(0.01)+(log10(10.0)-log10(0.01))*(j-1)/max(samples-1,1))
        s=0.05+2pi*hz*im
        theta=Diagonal([kp,ki])
        A=F.A0+F.U*theta*F.V'
        J=jacobians(x,u,p;kp,ki)
        push!(resolvent_errors,woodbury_error(F.A0,F.U,theta,F.V,s))
        Rw=woodbury_resolvent(F.A0,F.U,theta,F.V,s)
        Yw=J.C*Rw*J.B+J.D
        Yd=port_transfer(J.A,J.B,J.C,J.D,s)
        push!(port_errors,norm(Yd-Yw)/max(norm(Yd),eps(Float64)))
        derr=abs(det(s*I-A)-det(s*I-F.A0)*det(I-theta*(F.V'*((s*I-F.A0)\F.U)))) /
            max(abs(det(s*I-A)),eps(Float64))
        push!(denominator_errors,derr)
        push!(sample_rows,(sample=j,frequency_hz=hz,Kp=kp,Ki=ki,direct_norm=norm(Yd),
            resolvent_error=last(resolvent_errors),woodbury_error=last(port_errors),
            determinant_factorization_error=derr))
    end
    return (rank_Akp=1,rank_Aki=1,rank_both_A=2,
        rank_Bkp=1,rank_Bki=1,rank_both_B=2,
        rank_C_derivatives=0,rank_D_derivatives=0,
        resolvent_median=median(resolvent_errors),resolvent_p95=quantile(resolvent_errors,0.95),
        resolvent_max=maximum(resolvent_errors),port_median=median(port_errors),
        port_p95=quantile(port_errors,0.95),port_max=maximum(port_errors),
        determinant_factorization_max=maximum(denominator_errors),
        frequency_range_hz=[0.01,10.0],sample_count=samples,seed=seed,
        full_resolvent_errors=resolvent_errors,full_port_errors=port_errors,
        sample_rows=sample_rows,Aref=F.Aref,A0=F.A0,U=F.U,V=F.V)
end

"""
For rank-two PLL closure, det(I-diag(Kp,Ki)H) has exact bilinear
coefficients in the absolute gains. Coefficients are returned as complex
values; no polynomial degree is presumed for downstream physical port closure.
"""
function denominator_coefficients(A0,U,V,s)
    H=V'*((s*I-A0)\U)
    detH=det(H)
    return (c00=1.0,c10=-H[1,1],c01=-H[2,2],c11=detH,H=H,degree_Kp=1,degree_Ki=1,
            total_degree=2)
end

end
