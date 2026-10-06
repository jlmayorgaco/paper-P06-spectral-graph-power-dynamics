using LinearAlgebra

"Construct the exact common-angle gauge direction from frozen equilibria."
function gauge_vector(net, model)
    g=zeros(Float64,model.n_dynamic)
    for row in eachrow(model.state_map)
        b=Int(row.bus); firststate=Int(row.first)
        if row.kind=="SG"
            g[Int(row.last)]=1.0
        elseif row.kind=="GFL"
            op=net.gfl[b].op
            g[firststate+2]=1.0
            g[firststate+5]=op.x[7]
            g[firststate+6]=-op.x[6]
        else
            error("unknown dynamic component kind $(row.kind)")
        end
    end
    g
end

"Full finite spectrum plus an explicit global-angle gauge test."
function spectrum_with_gauge(net,rho,kp,ki;gauge_tol=1e-8)
    m=CollectiveModel.mixed_jacobian(net,rho,kp,ki)
    F=svd(m.Ared); g=gauge_vector(net,m)
    gauge_res=norm(m.Ared*g)/max(norm(m.Ared)*norm(g),eps())
    λ,R=eigen(m.Ared); ixg=argmin(abs.(λ)); gauge_pole=λ[ixg]
    # Deflate the exact common-angle direction before the finite eigensolve.
    # This avoids mixing the gauge with a nearby physical root.
    gauge_ok=gauge_res<=gauge_tol
    if gauge_ok
        qg=g/norm(g)
        Q=nullspace(transpose(qg))
        Abar=transpose(Q)*m.Ared*Q
        pλ,Rp=eigen(Abar)
        critical_idx=argmax(real.(pλ))
        right=Q*Rp[:,critical_idx]
        λL,Lp=eigen(adjoint(Abar))
        left=Q*Lp[:,argmin(abs.(λL .- conj(pλ[critical_idx])))]
    else
        Q=Matrix{Float64}(I,size(m.Ared,1),size(m.Ared,1))
        Abar=m.Ared
        pλ=λ
        critical_idx=argmax(real.(pλ))
        right=R[:,critical_idx]
        λL,Lp=eigen(adjoint(m.Ared))
        left=Lp[:,argmin(abs.(λL .- conj(pλ[critical_idx])))]
    end
    nonnormality=norm(left)*norm(right)/max(abs(dot(left,right)),eps())
    nearids=sortperm(abs.(pλ))[1:min(6,length(pλ))]
    threshold=max(1e-10,1e-9*maximum(F.S))
    nullity=count(x->x<=threshold,F.S)
    return (;model=m,lambda=λ,physical=pλ,spectral_abscissa=maximum(real.(pλ)),
        critical=pλ[critical_idx],critical_right=right,critical_left=left,
        critical_nonnormality=nonnormality,gauge_vector=g,gauge_residual=gauge_res,
        gauge_pole=gauge_ok ? 0.0+0im : gauge_pole,raw_nearest_pole=gauge_pole,
        gauge_detected=gauge_ok,singular_values=F.S,
        quotient_basis=Q,quotient_matrix=Abar,
        left_null=F.U[:,end],right_null=F.V[:,end],zero_nullity=nullity,
        nearzero=pλ[nearids],critical_index=critical_idx)
end

"Classify the global gauge and any distinct physical near-zero mode."
function zero_structure(net,rho,kp,ki;zero_tol=1e-8,near_tol=1e-4)
    s=spectrum_with_gauge(net,rho,kp,ki;gauge_tol=zero_tol)
    near=s.physical[argmin(abs.(s.physical))]
    alg=count(z->abs(z)<=zero_tol,s.physical)+(s.gauge_detected ? 1 : 0)
    A=s.model.Ared; g=s.gauge_vector; w=pinv(A)*g
    jordan_res=norm(A*w-g)/max(norm(g),eps())
    if !s.gauge_detected
        cls=abs(near)<=near_tol ? "NUMERICALLY_UNRESOLVED" : "NO_PHYSICAL_ZERO"
    elseif abs(near)<=zero_tol && alg>=2
        cls=s.zero_nullity>=2 ? "SEMISIMPLE_2D_ZERO_INCLUDING_GAUGE" : "DEFECTIVE_JORDAN"
    elseif abs(near)<=near_tol
        cls="SIMPLE_PHYSICAL_ZERO"
    else
        cls="NO_PHYSICAL_ZERO"
    end
    merge(s,(classification=cls,physical_nearzero=near,
        algebraic_multiplicity_estimate=alg,geometric_multiplicity=s.zero_nullity,
        generalized_vector_residual=jordan_res,generalized_vector=w,
        zero_tolerance=zero_tol,near_zero_tolerance=near_tol))
end

"Re-derive the stock isolated-PLL reference using each frozen terminal voltage."
function pll_reference(net,buses,kpmin,kpmax,kimin,kimax)
    rows=NamedTuple[]
    for (j,b) in enumerate(buses)
        p=net.gfl[b].op.parameters; V=abs(net.voltage[b]); τ=p.pll_tau
        r=1/(3τ); q=1/τ-2r
        kp=(2r-3τ*r^2)/V; ki=(r^2-2τ*r^3)/V
        A=[0.0 0.0 -V*ki/τ; 1.0 0.0 -V*kp/τ; 0.0 1.0 -1/τ]
        roots=eigvals(A)
        kpclip=clamp(kp,kpmin[j],kpmax[j]); kiclip=clamp(ki,kimin[j],kimax[j])
        push!(rows,(bus=b,V_pu=V,tau_s=τ,r_requested_s_inv=r,q_s_inv=q,
            kp_isolated=kp,ki_isolated=ki,kp_used_seed=kpclip,ki_used_seed=kiclip,
            kp_bound_hit=kpclip!=kp,ki_bound_hit=kiclip!=ki,
            max_seed_root_real_s_inv=maximum(real.(roots)),
            all_seed_roots_le_minus_r=maximum(real.(roots))<=-r+1e-7,
            max_decay_s_inv=1/(3τ),vieta_sum_s_inv=-1/τ,
            coefficient_identity_residual=maximum(abs.([τ,1.0,V*kp,V*ki] .-
                [τ,1.0,τ*(r^2+2r*q),τ*r^2*q])),
            vieta_bound_residual=abs(3r-1/τ)))
    end
    DataFrame(rows)
end

"Validate exact closure sensitivities against centered full-state eigenvalue differences."
function validate_sensitivities(net,rho,kp,ki,s0,dλ;step_rel=2e-5)
    n=length(rho); rows=NamedTuple[]
    function select_pole(rp,kpp,kip)
        λ=eigvals(CollectiveModel.mixed_jacobian(net,rp,kpp,kip).Ared)
        λ[argmin(abs.(λ.-s0))]
    end
    for group in (:rho,:kp,:ki), i in 1:n
        scale=group===:rho ? 2e-5 : step_rel*max(abs(group===:kp ? kp[i] : ki[i]),1.0)
        rp=copy(rho); rm=copy(rho); kpP=copy(kp); kpM=copy(kp); kiP=copy(ki); kiM=copy(ki)
        if group===:rho
            rp[i]+=scale; rm[i]-=scale
        elseif group===:kp
            kpP[i]+=scale; kpM[i]-=scale
        else
            kiP[i]+=scale; kiM[i]-=scale
        end
        fd=(select_pole(rp,kpP,kiP)-select_pole(rm,kpM,kiM))/(2scale)
        analytic=group===:rho ? dλ.ds_drho[i] : group===:kp ? dλ.ds_dKp[i] : dλ.ds_dKi[i]
        err=abs(analytic-fd)/max(abs(fd),1e-8)
        push!(rows,(parameter=String(group),index=i,analytic_real=real(analytic),
            analytic_imag=imag(analytic),finite_difference_real=real(fd),
            finite_difference_imag=imag(fd),relative_error=err,absolute_error=abs(analytic-fd),
            step=scale,closure_residual=dλ.closure_residual))
    end
    DataFrame(rows)
end

"Minimum normalized-energy gain correction with explicit box active-set updates."
function gain_correction(J,residual,x,lo,hi;scales=max.(abs.(x),1e-8))
    n=length(x); free=trues(n); dx=zeros(n); active=String[]; rem=copy(residual)
    for _ in 1:n+1
        idx=findall(free); isempty(idx) && break
        W=Diagonal(scales[idx].^2); A=J[:,idx]
        cand=W*transpose(A)*pinv(A*W*transpose(A))*rem
        viol=[k for k in eachindex(idx) if x[idx[k]]+cand[k]<lo[idx[k]] ||
              x[idx[k]]+cand[k]>hi[idx[k]]]
        if isempty(viol)
            dx[idx].=cand; break
        end
        score(k)=max((lo[idx[k]]-x[idx[k]]-cand[k])/scales[idx[k]],
                     (x[idx[k]]+cand[k]-hi[idx[k]])/scales[idx[k]],0.0)
        k=viol[argmax(score.(viol))]; j=idx[k]
        bound=cand[k]>0 ? hi[j] : lo[j]
        dx[j]=bound-x[j]; free[j]=false; push!(active,"$j:$bound")
        rem .-= J[:,j].*dx[j]
    end
    (;delta=dx,active_bounds=active,normalized_norm=norm(dx./scales),
        residual_after=J*dx-residual)
end

"Exhaustively enumerate vertices of the box LP with two collective inequalities."
function enumerate_surrogate(P,gamma,h,bS,bR;tol=1e-9)
    n=length(P); cands=NamedTuple[]
    for mask in 0:(2^n-1)
        upper=findall(i->((mask>>(i-1))&1)==1,1:n)
        base=zeros(n); base[upper].=1.0
        rs=bS-dot(gamma,base); rr=bR-dot(h,base); free=setdiff(1:n,upper)
        function accept(x,support,kind)
            all(x.>=-tol) && all(x.<=1+tol) || return
            dot(gamma,x)+tol>=bS && dot(h,x)+tol>=bR || return
            x=clamp.(x,0,1)
            push!(cands,(eps=x,cost=dot(P,x),support=support,case=kind,
                spectral_slack=dot(gamma,x)-bS,rocof_slack=dot(h,x)-bR,
                saturated=join(upper,":")))
        end
        accept(base,upper,"box_vertex")
        for i in free
            for (a,b,label) in ((gamma[i],rs,"spectral_edge"),(h[i],rr,"rocof_edge"))
                abs(a)>tol || continue
                x=copy(base); x[i]=b/a
                accept(x,sort(vcat(upper,i)),label)
            end
        end
        for ai in 1:length(free)-1, aj in ai+1:length(free)
            i,j=free[ai],free[aj]; D=gamma[i]*h[j]-gamma[j]*h[i]
            abs(D)>tol*max(abs(gamma[i]*h[j]),abs(gamma[j]*h[i]),1.0) || continue
            x=copy(base); x[i]=(rs*h[j]-rr*gamma[j])/D
            x[j]=(rr*gamma[i]-rs*h[i])/D
            accept(x,sort(vcat(upper,i,j)),"both_constraints_active")
        end
    end
    isempty(cands) && return (;feasible=false,eps=fill(NaN,n),cost=Inf,
        candidates=DataFrame(),status="INFEASIBLE_SURROGATE")
    sort!(cands,by=x->x.cost)
    table=DataFrame(bus_support=[join(c.support,":") for c in cands],
        case=[c.case for c in cands],retained_sg_MW=[c.cost for c in cands],
        spectral_slack=[c.spectral_slack for c in cands],rocof_slack=[c.rocof_slack for c in cands],
        upper_saturated=[c.saturated for c in cands])
    (;feasible=true,eps=cands[1].eps,cost=cands[1].cost,support=cands[1].support,
        candidates=table,status="SURROGATE_GLOBAL")
end

"Matrix determinant lemma radial update for the selected one/two-anchor direction."
function low_rank_boundary(net,epsbar,kp,ki,sigma_eff;omega=0.0)
    buses=sort(collect(keys(net.sg))); anchors=findall(>(1e-12),epsbar)
    sb=complex(-sigma_eff,omega)
    base=ClosureSpectrum.port_closure(net,sb,ones(length(epsbar)),kp,ki)
    C0=base.C; G=base.G
    isempty(anchors) && return (;status="NO_ANCHORS",roots=Float64[],
        determinant_residual=abs(det(C0)))
    cols=Int[]; Dblocks=Matrix{ComplexF64}[]
    for i in anchors
        b=buses[i]; pi=(2i-1):(2i)
        sg=CollectiveModel.port_admittance(net.sg[b].J,sb)
        op=net.gfl[b].op
        J=CollectiveModel.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=kp[i],ki=ki[i])
        gf=CollectiveModel.port_admittance(J,sb)
        append!(cols,pi); push!(Dblocks,epsbar[i].*(sg-gf))
    end
    U=zeros(ComplexF64,size(C0,1),length(cols)); W=zeros(ComplexF64,length(cols),size(C0,2))
    for (k,col) in enumerate(cols); U[col,k]=1; end
    for (a,i) in enumerate(anchors); W[2a-1:2a,:].=G[2i-1:2i,:]; end
    D=zeros(ComplexF64,length(cols),length(cols))
    for a in eachindex(anchors); D[2a-1:2a,2a-1:2a].=Dblocks[a]; end
    M=D*W*(C0\U); μ=eigvals(M)
    roots=[real(-inv(z)) for z in μ if abs(z)>1e-14 &&
        abs(imag(-inv(z)))<1e-8 && real(-inv(z))>=-1e-8]
    lemma=det(C0)*det(I+M); full=det(C0+U*D*W)
    (;status=isempty(roots) ? "NO_PHYSICAL_RADIAL_ROOT" : "RADIAL_ROOTS",
        roots=sort(roots),eigenvalues=μ,small_matrix=M,determinant_residual=abs(det(C0)),
        anchors=anchors,boundary_pole=sb,determinant_lemma_residual=abs(full-lemma))
end

