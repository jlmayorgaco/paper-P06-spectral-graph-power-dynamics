module PLLLowRank

using LinearAlgebra
using ...CollectiveModel: AnalyticGFLPLL

export pll_factors, woodbury_port, factor_audit

"The two gain rows share one phase-error covector in both A and B."
function pll_factors(op)
    p=op.parameters
    d=AnalyticGFLPLL.gain_derivatives(op.x,op.u,p)
    J0=AnalyticGFLPLL.jacobians(op.x,op.u,p;kp=0.0,ki=0.0)
    ex=Float64.(d.error_gradient_x)
    eu=Float64.(d.error_gradient_u)
    hp=zeros(9);hp[4]=1/p.pll_tau
    hi=zeros(9);hi[5]=1
    return (A0=J0.A,B0=J0.B,C=J0.C,D=J0.D,
        hp=hp,hi=hi,ex=ex,eu=eu,
        A_kp=d.A_kp,A_ki=d.A_ki,B_kp=d.B_kp,B_ki=d.B_ki)
end

"Exact one-scalar Woodbury GFL port, with independent Kp and Ki."
function woodbury_port(f,s,Kp,Ki)
    h=Kp.*f.hp+Ki.*f.hi
    R0=(s*I(9)-f.A0)\I(9)
    Rh=R0*h
    den=1-dot(f.ex,Rh)
    abs(den)>1e-14 || throw(SingularException(1))
    X=R0*f.B0+Rh*((transpose(f.ex)*R0*f.B0+transpose(f.eu))./den)
    f.D+f.C*X
end

function factor_audit(f)
    Ap=f.hp*transpose(f.ex);Ai=f.hi*transpose(f.ex)
    Bp=f.hp*transpose(f.eu);Bi=f.hi*transpose(f.eu)
    return (rank_Akp=rank(f.A_kp),rank_Aki=rank(f.A_ki),
        joint_rank=rank(hcat(f.A_kp,f.A_ki)),
        A_factor_error=max(norm(Ap-f.A_kp),norm(Ai-f.A_ki)),
        B_factor_error=max(norm(Bp-f.B_kp),norm(Bi-f.B_ki)))
end

end
