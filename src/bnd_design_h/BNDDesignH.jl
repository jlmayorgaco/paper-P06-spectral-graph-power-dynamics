module BNDDesignH

using LinearAlgebra
using Statistics
using ..BNDDesignG

const CM = BNDDesignG.CollectiveModel
const CS = BNDDesignG.ClosureSpectrum
const PORT_G_CACHE=Dict{String,Matrix{Float64}}()
const PLL_FACTOR_CACHE=Dict{Tuple{String,Int},Any}()

export quotient_model, jordan_audit, closure_coefficients, authority_vector,
       authority_derivative, single_authority, woodbury_audit, active_gain_basis,
       exact_self_energy_split, graph_alignment, direct_resolvent_peak,
       descriptor_finite_spectrum, nearest_physical_root

"Exact linear quotient by the analytically constructed uniform-angle gauge."
function quotient_model(net, epsv, kp, ki; gauge_tol=1e-11)
    rho = 1 .- Float64.(epsv)
    m = CM.mixed_jacobian(net, rho, kp, ki)
    g = BNDDesignG.gauge_vector(net, m)
    ng = norm(g)
    ng > 0 || error("empty rotational gauge vector")
    qg = g / ng
    gres = norm(m.Ared * qg) / max(norm(m.Ared), eps(Float64))
    Q = nullspace(reshape(qg, 1, :))
    Aq = transpose(Q) * m.Ared * Q
    gres <= gauge_tol || error("rotational gauge residual $gres exceeds $gauge_tol")
    (;model=m, gauge=g, gauge_residual=gres, Q, Aq, values=eigvals(Aq),
      finite_dimension=size(Aq,1), removed_dimension=size(m.Ared,1)-size(Aq,1))
end

function nearest_physical_root(net, epsv, kp, ki)
    q = quotient_model(net, epsv, kp, ki)
    j = argmin(abs.(q.values))
    (;lambda=q.values[j], index=j, quotient=q)
end

"Gauge and Jordan-chain diagnostics; multiplicity thresholds are reported, not hidden."
function jordan_audit(net, kp, ki; zero_tol=1e-7)
    q = quotient_model(net, zeros(length(kp)), kp, ki)
    A, Aq, g = q.model.Ared, q.Aq, q.gauge
    svq = svdvals(Aq); sva = svdvals(A)
    # Absolute s^-1 tolerance is reported and chosen below the observed
    # separation to the next singular direction; scaling by ||A|| would
    # misclassify several small-but-nonzero singular values as null.
    tq = 1e-8
    ta = 1e-8
    nullq = count(<=(tq), svq); nulla = count(<=(ta), sva)
    vals = eigvals(Aq)
    algq = count(z -> abs(z) <= zero_tol, vals)
    fullvals = eigvals(A)
    algfull = count(z -> abs(z) <= zero_tol, fullvals)
    Fq = svd(Aq)
    r = Fq.V[:,end]; l = Fq.U[:,end]
    # The least-squares generalized-vector equation tests whether the quotient
    # zero itself starts another Jordan chain.
    wq = pinv(Aq; rtol=1e-12) * r
    qchain = norm(Aq*wq-r) / max(norm(r), eps(Float64))
    wg = pinv(A; rtol=1e-12) * g
    gchain = norm(A*wg-g) / max(norm(g), eps(Float64))
    chain_length = gchain <= 1e-6 ? 2 : 1
    algfull_est = max(algfull, nulla + (chain_length > 1 ? 1 : 0))
    cls = nullq == 0 ? "NO_PHYSICAL_ZERO" :
          algq == 1 && nullq == 1 ? "SIMPLE_ZERO_AFTER_QUOTIENT" :
          algq > nullq ? "DEFECTIVE_AFTER_QUOTIENT" :
          algq == nullq && nullq > 1 ? "SEMISIMPLE_MULTIPLE_ZERO" : "UNRESOLVED"
    (;classification=cls, quotient_nullity=nullq, quotient_algebraic_multiplicity=algq,
      full_nullity=nulla, full_algebraic_multiplicity=algfull_est,
      raw_full_eigenvalue_count_within_tolerance=algfull,jordan_chain_length=chain_length,
      quotient_smallest_singular=minimum(svq), full_smallest_singular=minimum(sva),
      gauge_residual=q.gauge_residual, quotient_right_null=r, quotient_left_null=l,
      quotient_chain_residual=qchain, gauge_generalized_vector=wg,
      gauge_chain_residual=gchain, quotient=q)
end

"Taylor coefficients of det(C(s,epsilon)) about s=0 for one retained SG port."
function closure_coefficients(net, buses, i, kp, ki, e; h=1e-3)
    n=length(buses); 1<=i<=n || throw(BoundsError(buses,i))
    ep=zeros(n); ep[i]=e
    rho=1 .- ep
    dp=det(CS.port_closure(net, complex(h), rho, kp, ki).C)
    dm=det(CS.port_closure(net, complex(-h), rho, kp, ki).C)
    a1=(dp-dm)/(2h)
    a2=(dp+dm)/(2h^2)
    (;a1,a2,dp,dm,h,e)
end

function port_network_transfer(net,buses)
    get!(PORT_G_CACHE,net.root) do
        P=zeros(Float64,size(net.y_static,1),2length(buses))
        for (j,b) in enumerate(buses);P[2b-1:2b,2j-1:2j].=I(2);end
        transpose(P)*(net.y_static\P)
    end
end

function pll_factors(net,b)
    get!(PLL_FACTOR_CACHE,(net.root,b)) do
        op=net.gfl[b].op;p=op.parameters
        d=CM.AnalyticGFLPLL.gain_derivatives(op.x,op.u,p)
        J0=CM.AnalyticGFLPLL.jacobians(op.x,op.u,p;kp=0.0,ki=0.0)
        hp=zeros(9);hp[4]=1/p.pll_tau;hi=zeros(9);hi[5]=1
        (A0=J0.A,B0=J0.B,C=J0.C,D=J0.D,hp=hp,hi=hi,
         ex=Float64.(d.error_gradient_x),eu=Float64.(d.error_gradient_u))
    end
end

function woodbury_gfl_port(f,s,kp,ki)
    h=kp.*f.hp+ki.*f.hi
    R0=(s*I(9)-f.A0)\I(9);Rh=R0*h
    den=1-transpose(f.ex)*Rh
    X=R0*f.B0+Rh*((transpose(f.ex)*R0*f.B0+transpose(f.eu))/den)
    f.D+f.C*X
end

"The exact 20-port closure using the shared-scalar PLL Woodbury form."
function closure_with_gains(net,buses,s,rho,kp,ki)
    n=length(buses);G=port_network_transfer(net,buses)
    Y=zeros(ComplexF64,2n,2n)
    for (j,b) in enumerate(buses)
        pi=2j-1:2j
        Ys=CM.port_admittance(net.sg[b].J,s)
        Yg=woodbury_gfl_port(pll_factors(net,b),s,kp[j],ki[j])
        Y[pi,pi].=(1-rho[j]).*Ys+rho[j].*Yg
    end
    C=Matrix{ComplexF64}(I,2n,2n)+Y*G
    (;C,G,Y)
end

"Leading physical-root authority after factoring the exact common-angle zero."
function authority_vector(net, buses, kp, ki; e=1e-4, h=1e-3)
    [single_authority(net,buses,i,kp,ki;e=e,h=h) for i in eachindex(buses)]
end

function single_authority(net,buses,i,kp,ki;e=1e-4,h=1e-3)
    n=length(buses);rho=ones(n);pi=2i-1:2i
    a1eps=0.0+0.0im;a2=0.0+0.0im
    for (k,s) in enumerate((complex(h),complex(-h)))
        cl=closure_with_gains(net,buses,s,rho,kp,ki);C=cl.C
        b=buses[i];sg=CM.port_admittance(net.sg[b].J,s)
        gf=woodbury_gfl_port(pll_factors(net,b),s,kp[i],ki[i])
        dY=zeros(ComplexF64,2n,2n);dY[pi,pi].=sg-gf
        dC=dY*cl.G
        ddet=det(C)*tr(C\dC)
        if k==1
            a1eps+=ddet/(2h);a2+=det(C)/(2h^2)
        else
            a1eps-=ddet/(2h);a2+=det(C)/(2h^2)
        end
    end
    z=a1eps/a2
    (eltype(kp)<:Complex||eltype(ki)<:Complex) ? z : real(z)
end

"Symmetric gain derivatives in normalized log-box coordinates."
function authority_derivative(net, buses, kp, ki; relstep=1e-4, authority_h=1e-3)
    n=length(buses); z0=zeros(2n); rho=ones(n)
    base=authority_vector(net,buses,kp,ki;h=authority_h)
    J=zeros(n,2n); Jcs=similar(J);Jcoarse=similar(J);Jfine=similar(J)
    for i in 1:n
        eps_plus=0.0+0.0im;eps_minus=0.0+0.0im
        a2=0.0+0.0im;d_a2=zeros(ComplexF64,2n)
        d_a1eps=zeros(ComplexF64,2n)
        for (q,s) in enumerate((complex(authority_h),complex(-authority_h)))
            cl=CS.port_closure(net,s,rho,kp,ki;derivatives=true);C=cl.C
            Ceps=-cl.d_rho[i]
            ddet=det(C);R_Ceps=C\Ceps
            de=ddet*tr(R_Ceps)
            if q==1;eps_plus=de;a2+=ddet/(2authority_h^2)
            else;eps_minus=de;a2+=ddet/(2authority_h^2) end
            Dparts=vcat(cl.d_kp,cl.d_ki)
            for j in 1:2n
                Cp=Dparts[j]
                R_Cp=C\Cp
                dp=ddet*tr(R_Cp)
                d_a2[j]+=dp/(2authority_h^2)
                # d^2(det C)/(d epsilon_i d K_j), including the direct
                # mixed port derivative only when both parameters share bus i.
                Cep=(j==i || j==n+i) ? -Cp : zeros(ComplexF64,size(C))
                mixed=ddet*(tr(R_Ceps)*tr(R_Cp)+tr(C\Cep)-tr(C\(Cp*R_Ceps)))
                d_a1eps[j]+= (q==1 ? 1 : -1)*mixed/(2authority_h)
            end
        end
        a1eps=(eps_plus-eps_minus)/(2authority_h)
        Ai=real(a1eps/a2)
        for j in 1:2n
            J[i,j]=real((d_a1eps[j]-Ai*d_a2[j])/a2)
        end
    end
    # Independent symmetric finite differences validate the analytic mixed
    # derivative. The two step sizes also expose cancellation near zero.
    for j in 1:2n
        x0=j<=n ? kp[j] : ki[j-n]
        cs=1e-20*max(abs(x0),1.0)
        ppc=ComplexF64.(kp);ipc=ComplexF64.(ki)
        if j<=n;ppc[j]+=im*cs else;ipc[j-n]+=im*cs end
        Jcs[:,j].=imag.(authority_vector(net,buses,ppc,ipc;h=authority_h))./cs
        for (factor,dest) in ((relstep,Jcoarse),(relstep/2,Jfine))
            step=factor*x0;pp=copy(kp);pm=copy(kp);ip=copy(ki);imv=copy(ki)
            if j<=n;pp[j]+=step;pm[j]-=step else;ip[j-n]+=step;imv[j-n]-=step end
            dest[:,j].=(authority_vector(net,buses,pp,ip;h=authority_h)-
                        authority_vector(net,buses,pm,imv;h=authority_h))./(2step)
        end
    end
    rel=abs.(J-Jcs)./max.(abs.(J),1e-8)
    finite=filter(isfinite,vec(rel))
    (;base,J,Jcs,Jcoarse,Jfine,relative_error=rel,
      median_relative_error=isempty(finite) ? NaN : median(finite),
      p95_relative_error=isempty(finite) ? NaN : quantile(finite,0.95),
      normalized_jacobian=J*Diagonal(vcat(kp,ki)), z0)
end

function woodbury_audit(net,buses)
    rows=NamedTuple[]
    for b in buses
        op=net.gfl[b].op; p=op.parameters
        d=CM.AnalyticGFLPLL.gain_derivatives(op.x,op.u,p)
        ex=Float64.(d.error_gradient_x); eu=Float64.(d.error_gradient_u)
        hp=zeros(9);hp[4]=1/p.pll_tau; hi=zeros(9);hi[5]=1
        epA=hp*transpose(ex); eiA=hi*transpose(ex)
        epB=hp*transpose(eu); eiB=hi*transpose(eu)
        aerr=max(norm(epA-d.A_kp),norm(eiA-d.A_ki))
        berr=max(norm(epB-d.B_kp),norm(eiB-d.B_ki))
        push!(rows,(bus=b,rank_Kp=rank(d.A_kp),rank_Ki=rank(d.A_ki),
          joint_rank=rank(hcat(d.A_kp,d.A_ki)),A_factor_error=aerr,
          B_factor_error=berr,
          representation="EXACT_ONE_SCALAR_WOODBURY"))
    end
    rows
end

function active_gain_basis(J; energy=0.999)
    F=svd(Matrix{Float64}(J))
    den=sum(abs2,F.S); den>0 || return (;V=zeros(size(J,2),0),singular_values=F.S,
        rank=0,captured=0.0,energy)
    cumulative=cumsum(abs2.(F.S))./den
    r=findfirst(>=(energy),cumulative)
    r===nothing && (r=length(F.S))
    (;V=F.V[:,1:r],singular_values=F.S,rank=r,captured=cumulative[r],energy)
end

"Exact determinant-derivative split into nodal and off-diagonal network pathways."
function exact_self_energy_split(net,buses,kp,ki;h=1e-3)
    n=length(buses); rows=NamedTuple[]; pathrows=NamedTuple[]
    ep=zeros(n); rho=ones(n)
    for i in 1:n
        dD_parts=zeros(Float64,2,2)
        for (q,s) in enumerate((complex(h),complex(-h)))
            cl=CS.port_closure(net,s,rho,kp,ki)
            C=cl.C; G=cl.G
            Gd=zeros(ComplexF64,size(G))
            for j in 1:n
                jj=2j-1:2j; Gd[jj,jj].=G[jj,jj]
            end
            b=buses[i]; pi=2i-1:2i
            sg=CM.port_admittance(net.sg[b].J,s)
            op=net.gfl[b].op
            Jg=CM.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=kp[i],ki=ki[i])
            gf=CM.port_admittance(Jg,s)
            dY=zeros(ComplexF64,2n,2n); dY[pi,pi].=sg-gf
            for (part,Gpart,col) in ((1,Gd,1),(2,G-Gd,2))
                dC=dY*Gpart
                ddet=det(C)*tr(C\dC)
                dD_parts[q,col]=real(ddet)
            end
        end
        a1direct=(dD_parts[1,1]-dD_parts[2,1])/(2h)
        a1self=(dD_parts[1,2]-dD_parts[2,2])/(2h)
        cp=closure_coefficients(net,buses,i,kp,ki,1e-5;h=h)
        a2=cp.a2
        A_direct=real(a1direct/a2); A_self=real(a1self/a2)
        total=A_direct+A_self
        push!(rows,(bus=buses[i],A_direct=A_direct,A_self_energy=A_self,
          A_total=total,self_energy_share=abs(A_self)/max(abs(total),eps()),
          split_status="EXACT_ADJUGATE_DERIVATIVE_G_DIAG_OFFDIAG"))
        # Exact port-space Feshbach partition for anchor k. Its Schur self-
        # energy is Gamma=-T_kr*T_rr^{-1}*T_rk; each (l,m) term is retained.
        cl0=CS.port_closure(net,complex(h),rho,kp,ki);C=cl0.C
        kiport=2i-1:2i; rest=setdiff(1:2n,collect(kiport))
        Tkr=C[kiport,rest];Trk=C[rest,kiport];Trr=C[rest,rest]
        X=Trr\Trk;Gamma=-(Tkr*X)
        Tinv=Trr\I(size(Trr,1)); kappa=0.0
        for l in 1:n, m in 1:n
            (l==i||m==i)&&continue
            pl=2l-1:2l;pm=2m-1:2m
            rl=[findfirst(==(z),rest) for z in pl];rm=[findfirst(==(z),rest) for z in pm]
            gl=-C[kiport,pl]*(Tinv[rl,rm])*C[pm,kiport]
            val=norm(gl);kappa+=val
            push!(pathrows,(anchor_bus=buses[i],via_bus_l=buses[l],via_bus_m=buses[m],
              pathway_norm=val,pathway_trace_real=real(tr(gl)),
              Gamma_frobenius=norm(Gamma),pathway_cancellation_ratio=NaN))
        end
        for rowidx in length(pathrows)-((n-1)^2)+1:length(pathrows)
            pathrows[rowidx]=merge(pathrows[rowidx],(pathway_cancellation_ratio=kappa/max(norm(Gamma),eps()),))
        end
    end
    (;authority=rows,pathways=pathrows)
end

function graph_alignment(net,buses,dispatch,authority,Vactive,kp,ki)
    gf=BNDDesignG.graph_features(net,buses,dispatch,authority)
    n=length(buses); Bg=zeros(Float64,2n,2gf.rank)
    gf.rank>0 && (Bg[1:n,1:gf.rank].=gf.basis; Bg[n+1:2n,gf.rank+1:2gf.rank].=gf.basis)
    Bg=Matrix(qr(Bg).Q)[:,1:rank(Bg)]
    singular=svdvals(transpose(Vactive)*Bg)
    angles=acos.(clamp.(singular,0.0,1.0))
    P=Bg*transpose(Bg)
    J=Vactive
    captured=norm(J*transpose(J)*P)^2/max(norm(J*transpose(J))^2,eps())
    (;graph=gf,Bg,principal_angles_rad=angles,captured_energy=captured)
end

function direct_resolvent_peak(A,sigma;scale=1.0)
    lambda=eigvals(A);alpha=maximum(real.(lambda))
    if alpha>=-sigma-1e-11
        return (;status="NOMINAL_NOT_LEFT_OF_BOUNDARY",peak=Inf,omega_peak=NaN,
          beta_star=0.0,frequency=Float64[],sigma_min=Float64[],refinements=0)
    end
    As=A+sigma*I
    wmax=max(1.0,4maximum(abs.(imag.(lambda)),init=0.0),20maximum(abs.(lambda),init=0.0))
    grid=unique(vcat(0.0,10.0.^range(-4,log10(wmax),length=72),collect(range(0.0,wmax,length=88))))
    f(w)=scale*opnorm((im*w*I-As)\I,2)
    vals=f.(grid);ix=argmax(vals);peak=vals[ix];wp=grid[ix];refinements=0
    for i in 2:length(grid)-1
        vals[i]>=vals[i-1]&&vals[i]>=vals[i+1]||continue
        a=grid[i-1];b=grid[i+1]
        c=b-(b-a)*0.6180339887498949;d=a+(b-a)*0.6180339887498949
        fc=f(c);fd=f(d)
        for _ in 1:32
            if fc>fd
                b=d;d=c;fd=fc;c=b-(b-a)*0.6180339887498949;fc=f(c)
            else
                a=c;c=d;fc=fd;d=a+(b-a)*0.6180339887498949;fd=f(d)
            end
        end
        w=(a+b)/2;v=f(w);refinements+=1
        if v>peak;peak=v;wp=w;end
    end
    (;status="EVALUATED",peak,omega_peak=wp,beta_star=inv(peak),frequency=grid,
      sigma_min=1.0./max.(vals,eps(Float64)),refinements)
end

"Finite generalized eigenvalues of the unreduced state/algebraic pencil."
function descriptor_finite_spectrum(m)
    n=size(m.A,1); p=size(m.Gy,1)
    Ad=zeros(Float64,n+p,n+p); Ed=zeros(Float64,n+p,n+p)
    # E xdot = A x + B y, 0 = C x + Gy y. Writing this as
    # (A_d - s E_d)[x;y]=0 gives the exact Schur-reduced Ared.
    Ad[1:n,1:n].=m.A; Ad[1:n,n+1:end].=m.B
    Ad[n+1:end,1:n].=-m.C; Ad[n+1:end,n+1:end].=-m.Gy
    Ed[1:n,1:n].=I(n)
    z=eigvals(Ad,Ed)
    filter(x->isfinite(real(x))&&isfinite(imag(x)),z)
end

end
