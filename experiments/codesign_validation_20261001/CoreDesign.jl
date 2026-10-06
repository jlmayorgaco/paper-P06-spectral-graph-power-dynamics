module CoreDesign

using LinearAlgebra, Statistics, DataFrames, TOML
import ..ExpP
import ..LinearSecurity
import ..FiniteWindow

const N=ExpP.PDExactDesignN
const SIGMA=0.05
const BETA_REQ=1.6991206999182038e-6
const FREQUENCY_LIMIT=0.5
const ROCOF_LIMIT=0.5

export decode, encode, design_eval, robust_beta_sampled, solve_fixed_support,
       constraint_jacobian, kkt_audit

function decode(ctx,support,y)
    ns=length(support); length(y)==ns+20 || throw(DimensionMismatch())
    epsv=Float64.(y[1:ns]); kp=ctx.kpmin .+ Float64.(y[ns+1:ns+10]).*(ctx.kpmax.-ctx.kpmin)
    ki=ctx.kimin .+ Float64.(y[ns+11:ns+20]).*(ctx.kimax.-ctx.kimin)
    epsv,kp,ki
end

function encode(ctx,epsv,kp,ki)
    vcat(Float64.(epsv),(kp.-ctx.kpmin)./(ctx.kpmax.-ctx.kpmin),
         (ki.-ctx.kimin)./(ctx.kimax.-ctx.kimin))
end

"Sampled ExpG full-block beta radius for the exact frozen shift.

The optimization evaluator uses a deterministic adaptive candidate set made
from zero, a low-frequency mesh, and lightly damped pole frequencies. This is
an observed estimate. The final candidate is rechecked with ExpG's full
resolvent refinement and interval lower bound before any robust claim.
"
function robust_beta_sampled(A,lambda;sigma=SIGMA,scale=1.0)
    maximum(real.(lambda)) < -sigma || return (beta=0.0,omega=NaN,
        peak=Inf,samples=0,status="NOMINAL_NOT_LEFT_OF_SHIFT")
    As=ComplexF64.(A+sigma*I)
    # The local active-set corrector monitors zero and the frequencies of the
    # three least-damped modes. Full ExpG frequency refinement is mandatory
    # for the frozen result and is never inferred from this small set.
    slow=sortperm(real.(lambda);rev=true)[1:min(3,length(lambda))]
    grid=unique(vcat(0.0,abs.(imag.(lambda[slow]))))
    best=Inf;wb=NaN
    for w in grid
        σ=minimum(svdvals(im*Float64(w)*I(size(As,1))-As))/scale
        if σ<best;best=σ;wb=Float64(w);end
    end
    (;beta=best,omega=wb,peak=1/best,samples=length(grid),status="SAMPLED_ESTIMATE")
end

function _freq_outputs(ctx,m,rho)
    yload=LinearSecurity.load_input_vector(ctx,m,16;system_base_mva=100.0)
    Bfull=-m.B*(m.Gy\yload)
    g=N.gauge_vector(m);Q=nullspace(reshape(g/norm(g),1,:))
    A=transpose(Q)*m.Ared*Q;B=transpose(Q)*Bfull
    Cdot=-(m.Gy\(m.C*m.Ared));Ddot=-(m.Gy\(m.C*Bfull))
    Cbus=zeros(10,size(m.Ared,1));Dbus=zeros(10)
    for bus in 30:39
        k=bus-29;v=ctx.net.voltage[bus];row=zeros(1,size(m.Gy,1))
        row[1,2bus-1]=-imag(v)/(abs2(v)*2pi);row[1,2bus]=real(v)/(abs2(v)*2pi)
        Cbus[k,:].=(row*Cdot)[:];Dbus[k]=(row*Ddot)[1]
    end
    (;A,B,C_bus=Cbus*Q,D_bus=Dbus,Q,gauge_residual=norm(m.Ared*g)/max(norm(m.Ared)*norm(g),eps()))
end

"Evaluate the complete physical constraints at a fixed support and design."
function design_eval(ctx,support,epsv,kp,ki;disturbance_MW=100.0,
                     rocof_window=nothing,dt_s=0.05,horizon_s=30.0)
    rho=ones(10)
    for (j,b) in enumerate(support);rho[b-29]=1-Float64(epsv[j]);end
    sp=N.spectrum(ctx,rho,kp,ki); alpha=sp.alpha
    beta=robust_beta_sampled(transpose(sp.quotient)*sp.model.Ared*sp.quotient,
                             sp.lambda)
    tm=FiniteWindow.design_metrics(ctx,sp.model,rho;load_bus=16,
        disturbance_MW,windows=(0.5,0.2,1.0,2.0),dt_s,horizon_s,
        gauge_vector=N.gauge_vector)
    fpeak=only(filter(r->r.window_s==0.5,tm.metrics))
    finf=maximum(abs.(tm.signals.F_inf_Hz))
    vals=Float64[(alpha+SIGMA)/0.01,
        (BETA_REQ-beta.beta)/BETA_REQ,
        (finf-FREQUENCY_LIMIT)/FREQUENCY_LIMIT,
        (fpeak.F_peak_Hz-FREQUENCY_LIMIT)/FREQUENCY_LIMIT]
    if rocof_window!==nothing
        row=only(filter(r->r.window_s==Float64(rocof_window),tm.metrics))
        push!(vals,(row.R_peak_Hz_s-ROCOF_LIMIT)/ROCOF_LIMIT)
    end
    (;g=vals,alpha,beta=beta.beta,beta_observed=beta,Finf=finf,
      Fpeak=fpeak.F_peak_Hz,Rwindows=tm.metrics,signals=tm.signals,
      lambda=sp.lambda,Q=sp.quotient,
      quotient=transpose(sp.quotient)*sp.model.Ared*sp.quotient,
      model=sp.model,rho,epsilon=Float64.(epsv),Kp=Float64.(kp),Ki=Float64.(ki),
      frequency_condition=tm.signals.condition_A,
      eigen_condition=tm.signals.eigenvector_condition)
end

function _bounds(ctx,support)
    ns=length(support);lo=vcat(fill(1e-5,ns),zeros(20));hi=vcat(fill(1-1e-5,ns),ones(20))
    lo,hi
end

function _robust_gradient(ctx,support,y,ev;step=2e-4)
    ω=Float64(ev.beta_observed.omega);Aq=ev.quotient;Q=ev.Q
    M=ComplexF64.(im*ω*I(size(Aq,1))-(Aq+SIGMA*I))
    F=svd(M);u=F.U[:,end];v=F.V[:,end]
    gap=length(F.S)>1 ? F.S[end-1]-F.S[end] : Inf
    ns=length(y);db=zeros(ns)
    for col in 1:ns
        h=step;yp=copy(y);ym=copy(y)
        yp[col]+=h;ym[col]-=h
        lo,hi=_bounds(ctx,support)
        plus=yp[col]<=hi[col];minus=ym[col]>=lo[col]
        if plus && minus
            ep,kp,ki=decode(ctx,support,yp);rp=ones(10)
            for (j,b) in enumerate(support);rp[b-29]=1-ep[j];end
            mp=N.descriptor(ctx,rp,kp,ki);Ap=transpose(Q)*mp.Ared*Q
            ep,kp,ki=decode(ctx,support,ym);rm=ones(10)
            for (j,b) in enumerate(support);rm[b-29]=1-ep[j];end
            mm=N.descriptor(ctx,rm,kp,ki);Am=transpose(Q)*mm.Ared*Q
            Az=(Ap-Am)/(2h)
        elseif plus
            ep,kp,ki=decode(ctx,support,yp);rp=ones(10)
            for (j,b) in enumerate(support);rp[b-29]=1-ep[j];end
            mp=N.descriptor(ctx,rp,kp,ki);Ap=transpose(Q)*mp.Ared*Q
            Az=(Ap-Aq)/h
        elseif minus
            ep,kp,ki=decode(ctx,support,ym);rm=ones(10)
            for (j,b) in enumerate(support);rm[b-29]=1-ep[j];end
            mm=N.descriptor(ctx,rm,kp,ki);Am=transpose(Q)*mm.Ared*Q
            Az=(Aq-Am)/h
        else
            continue
        end
        # Envelope derivative of sigma_min(iωI-Aσ); beta_star is that
        # minimum singular value for the frozen full-state block.
        db[col]=real(dot(u,-Az*v))
    end
    (;db,gap,omega=ω,singular_gap=gap)
end

function _dc_frequency(ctx,support,y,d)
    epsv,kp,ki=decode(ctx,support,y);rho=ones(10)
    for (j,b) in enumerate(support);rho[b-29]=1-epsv[j];end
    m=N.descriptor(ctx,rho,kp,ki);red=_freq_outputs(ctx,m,rho)
    h=real.(red.C_bus*(-(red.A\red.B)).+red.D_bus).*Float64(d)
    maximum(abs.(h))
end

function _window_measurement(ctx,support,y,d;windows=(0.5,),dt_s=0.05,horizon_s=30.0)
    epsv,kp,ki=decode(ctx,support,y);rho=ones(10)
    for (j,b) in enumerate(support);rho[b-29]=1-epsv[j];end
    m=N.descriptor(ctx,rho,kp,ki)
    FiniteWindow.design_metrics(ctx,m,rho;load_bus=16,disturbance_MW=d,
        windows,dt_s,horizon_s,gauge_vector=N.gauge_vector)
end

"Finite-difference only the active event/robust rows; modal row is analytic."
function constraint_jacobian(ctx,support,y,ev;rocof_window=nothing,
                              disturbance_MW=100.0,dt_s=0.05,
                              horizon_s=30.0,fd_step=5e-4,active_tol=0.5)
    n=length(y);m=length(ev.g);J=zeros(m,n)
    active=findall(ev.g .>= -active_tol)
    # Full-spectrum active pole sensitivity is analytic within this support.
    rho=copy(ev.rho)
    if 1 in active
        j=argmax(real.(ev.lambda))
        try
            ds=N.simple_mode_sensitivities(ctx,rho,ev.Kp,ev.Ki;mode=j)
            for (k,bus) in enumerate(support);J[1,k]=-real(ds.rho[bus-29])/0.01;end
            off=length(support)
            for i in 1:10
                J[1,off+i]=real(ds.Kp[i])*(ctx.kpmax[i]-ctx.kpmin[i])/0.01
                J[1,off+10+i]=real(ds.Ki[i])*(ctx.kimax[i]-ctx.kimin[i])/0.01
            end
        catch err
            rethrow(err) # a failed modal derivative must never become a zero row
        end
    end
    baseeps,basekp,baseki=decode(ctx,support,y)
    dc_done=false
    for col in 1:n
        # The exact DC gain is independent of PLL PI gains, so its derivative
        # is zero outside the retained-SG coordinates. This also avoids model
        # evaluations for inactive rows far from their constraints.
        needs=filter(r->!(r in (1,2,3,4,5)),active)
        # The peak-frequency constraint has the same gradient as the exact
        # DC output when its measured maximum is the steady plateau.
        if 3 in active && col<=length(support)
            h=fd_step;yp=copy(y);ym=copy(y);yp[col]+=h;ym[col]-=h
            lo,hi=_bounds(ctx,support);plus=yp[col]<=hi[col];minus=ym[col]>=lo[col]
            if plus&&minus
                fp=_dc_frequency(ctx,support,yp,disturbance_MW)
                fm=_dc_frequency(ctx,support,ym,disturbance_MW)
                J[3,col]=(fp-fm)/(2h*FREQUENCY_LIMIT)
            elseif plus
                fp=_dc_frequency(ctx,support,yp,disturbance_MW)
                J[3,col]=(fp-ev.Finf)/(h*FREQUENCY_LIMIT)
            elseif minus
                fm=_dc_frequency(ctx,support,ym,disturbance_MW)
                J[3,col]=(ev.Finf-fm)/(h*FREQUENCY_LIMIT)
            end
        end
        if 4 in active && abs(ev.Fpeak-ev.Finf)<=1e-7
            J[4,col]=J[3,col]
        end
        if (4 in active && abs(ev.Fpeak-ev.Finf)>1e-7) || (rocof_window!==nothing && 5 in active)
            h=fd_step;yp=copy(y);ym=copy(y);yp[col]+=h;ym[col]-=h
            lo,hi=_bounds(ctx,support);plus=yp[col]<=hi[col];minus=ym[col]>=lo[col]
            if (4 in active && abs(ev.Fpeak-ev.Finf)>1e-7)
                fp=plus ? only(_window_measurement(ctx,support,yp,disturbance_MW;
                    windows=(0.5,),dt_s,horizon_s).metrics).F_peak_Hz : ev.Fpeak
                fm=minus ? only(_window_measurement(ctx,support,ym,disturbance_MW;
                    windows=(0.5,),dt_s,horizon_s).metrics).F_peak_Hz : ev.Fpeak
                J[4,col]=(fp-fm)/((plus&&minus ? 2h : h)*FREQUENCY_LIMIT)
            end
            if rocof_window!==nothing && 5 in active
                fp=plus ? only(_window_measurement(ctx,support,yp,disturbance_MW;
                    windows=(Float64(rocof_window),),dt_s,horizon_s).metrics).R_peak_Hz_s :
                    only(filter(r->r.window_s==Float64(rocof_window),ev.Rwindows)).R_peak_Hz_s
                fm=minus ? only(_window_measurement(ctx,support,ym,disturbance_MW;
                    windows=(Float64(rocof_window),),dt_s,horizon_s).metrics).R_peak_Hz_s :
                    only(filter(r->r.window_s==Float64(rocof_window),ev.Rwindows)).R_peak_Hz_s
                J[5,col]=(fp-fm)/((plus&&minus ? 2h : h)*ROCOF_LIMIT)
            end
        end
    end
    if 2 in active && isfinite(ev.beta_observed.omega) && ev.beta>0
        bg=_robust_gradient(ctx,support,y,ev;step=fd_step)
        J[2,:].=-bg.db./BETA_REQ
    end
    # PLL gains have exactly zero DC output derivative for this frozen
    # SimpleGFLDC equilibrium; encode the proven identity directly.
    if 3 in active
        J[3,length(support)+1:end].=0.0
    end
    (;J,active_rows=active,analytic_modal=1 in active)
end

function _solve_qp(H,c,A,b;tol=1e-9,maxiter=250)
    n=length(c);d=zeros(n)
    minimum(b-A*d)>=-1e-8 || return (;d,status="LINEARIZED_INFEASIBLE",W=Int[],mu=Float64[],residual=Inf)
    W=findall(abs.(b-A*d).<=1e-10)
    for _ in 1:maxiter
        if isempty(W)
            p=-(H\(H*d+c));mu=Float64[]
        else
            Aw=A[W,:];K=[H transpose(Aw);Aw zeros(length(W),length(W))]
            sol=try K\vcat(-(H*d+c),zeros(length(W))) catch; pinv(K)*vcat(-(H*d+c),zeros(length(W))) end
            p=sol[1:n];mu=sol[n+1:end]
        end
        if norm(p,Inf)<=tol
            if isempty(mu)||minimum(mu)>=-1e-8
                stat=norm(H*d+c+(isempty(W) ? zeros(n) : transpose(A[W,:])*mu),Inf)
                return (;d,status="QP_OPTIMAL",W,mu,residual=stat)
            end
            deleteat!(W,argmin(mu));continue
        end
        α=1.0;block=0;inW=Set(W)
        for r in eachindex(b)
            r in inW && continue
            den=dot(view(A,r,:),p)
            if den>tol
                ar=(b[r]-dot(view(A,r,:),d))/den
                if ar<α;α=max(0.0,ar);block=r;end
            end
        end
        d .+= α.*p
        block!=0 && !(block in W) && push!(W,block)
    end
    (;d,status="QP_ITERATION_LIMIT",W,mu=Float64[],residual=Inf)
end

function _cost(ctx,support,epsv)
    sum(ctx.power[b-29]*epsv[j] for (j,b) in enumerate(support))
end

"Feasible-start active-set SQP on one fixed physical architecture.

The objective is linear, variables and trust region are normalized, every
accepted point is re-evaluated with the complete finite spectrum, sampled
ExpG resolvent radius, and finite-window phase response.
"
function solve_fixed_support(ctx,support,eps0,kp0,ki0;maxiter=12,
    trust_radius=0.05,rocof_window=nothing,disturbance_MW=100.0,
    dt_s=0.05,horizon_s=30.0,
    guard=2e-7,checkpoint_path=nothing,restoration=true,
    restoration_maxiter=4)
    s=sort!(unique(Int.(support))); y=encode(ctx,eps0,kp0,ki0)
    lo,hi=_bounds(ctx,s); all(lo .<= y .<= hi) || error("seed outside support bounds")
    ev=design_eval(ctx,s,eps0,kp0,ki0;disturbance_MW,rocof_window,dt_s,horizon_s)
    if maximum(ev.g)>guard
        return (;status="INFEASIBLE_SEED",support=s,y,epsilon=eps0,Kp=kp0,Ki=ki0,
            retained_SG_MW=_cost(ctx,s,eps0),evaluation=ev,history=DataFrame())
    end
    c=zeros(length(y)); for (j,b) in enumerate(s);c[j]=ctx.power[b-29]/1000;end
    H=Matrix{Float64}(I,length(y),length(y));history=NamedTuple[];accepted=0;term="MAX_ITERATIONS"
    for it in 1:maxiter
        jac=constraint_jacobian(ctx,s,y,ev;rocof_window,disturbance_MW,dt_s,horizon_s)
        ns=length(y);Id=Matrix{Float64}(I,ns,ns)
        # g(y)<=0, box and trust constraints. Keep a small interior buffer.
        A=vcat(jac.J,Id,-Id,Id,-Id)
        b=vcat(-ev.g,hi.-y,y.-lo,fill(trust_radius,length(y)),
            fill(trust_radius,length(y)))
        qp=_solve_qp(H,c,A,b)
        push!(history,(;iteration=it,cost=_cost(ctx,s,decode(ctx,s,y)[1]),
            max_constraint=maximum(ev.g),alpha=ev.alpha,beta=ev.beta,Finf=ev.Finf,
            Fpeak=ev.Fpeak,active_rows=join(jac.active_rows,";"),
            qp_status=qp.status,step_inf=norm(qp.d,Inf),accepted=false,
            line_search=0.0,stationarity_qp=qp.residual,trust_radius))
        println("SQP support=",join(s,":")," d=",disturbance_MW," it=",it,
            " J=",history[end].cost," maxg=",history[end].max_constraint,
            " alpha=",ev.alpha," beta=",ev.beta," trust=",trust_radius)
        if checkpoint_path!==nothing
            _checkpoint(checkpoint_path,s,y,ctx,ev,disturbance_MW,"ITERATION_START")
        end
        if qp.status!="QP_OPTIMAL";term=qp.status;break;end
        if norm(qp.d,Inf)<2e-6;term="QP_STATIONARY";break;end
        e0,k0p,k0i=decode(ctx,s,y);J0=_cost(ctx,s,e0);committed=false
        for a in (1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125)
            yt=y+a*qp.d
            all(lo .<= yt .<= hi) || continue
            et,kpt,kit=decode(ctx,s,yt)
            evt=design_eval(ctx,s,et,kpt,kit;disturbance_MW,
                rocof_window,dt_s,horizon_s)
            # Curved active boundaries (especially the small-gain radius near
            # beta_req) can invalidate a tangent SQP step even when its
            # first-order prediction is feasible. Correct such trials with
            # Newton projections onto the exact currently-near constraints;
            # every projected point is re-evaluated before acceptance.
            if restoration && maximum(evt.g)>guard
                for _ in 1:restoration_maxiter
                    maximum(evt.g)<=guard && break
                    rows=findall(evt.g .>= -0.05)
                    # The resolvent radius is defined as zero outside the
                    # nominally stable shifted domain. When a trial crosses
                    # the modal boundary, restore that boundary first; its
                    # modal derivative is well defined, while a beta
                    # derivative is not.
                    if evt.alpha >= -SIGMA && 2 in rows
                        deleteat!(rows,findfirst(==(2),rows))
                    end
                    isempty(rows) && break
                    jr=constraint_jacobian(ctx,s,yt,evt;rocof_window,
                        disturbance_MW,dt_s,horizon_s,fd_step=2e-4,
                        active_tol=0.05).J[rows,:]
                    residual=evt.g[rows].-(-guard/2)
                    correction=-(pinv(jr)*residual)
                    norm(correction,Inf)<=max(0.02,2trust_radius) || break
                    ycorr=yt+correction
                    all(lo .<= ycorr .<= hi) || break
                    yt=ycorr
                    et,kpt,kit=decode(ctx,s,yt)
                    evt=design_eval(ctx,s,et,kpt,kit;disturbance_MW,
                        rocof_window,dt_s,horizon_s)
                end
            end
            cost=_cost(ctx,s,et)
            if maximum(evt.g)<=1e-8 && cost<J0-1e-8
                y=yt;ev=evt;accepted+=1;committed=true
                history[end]=merge(history[end],(;cost,max_constraint=maximum(ev.g),
                    alpha=ev.alpha,beta=ev.beta,Finf=ev.Finf,Fpeak=ev.Fpeak,
                    accepted=true,line_search=a));break
            end
        end
        if !committed
            trust_radius*=0.5
            if trust_radius<1e-4;term="NO_FEASIBLE_DESCENT";break;end
        end
        checkpoint_path!==nothing && _checkpoint(checkpoint_path,s,y,ctx,ev,
            disturbance_MW,committed ? "ACCEPTED" : term)
    end
    e,kp,ki=decode(ctx,s,y)
    (;status=term,support=s,y,epsilon=e,rho=1 .- map(i->begin
        v=1.0; j=findfirst(==(i),s); j===nothing ? 0.0 : e[j]
    end,30:39),Kp=kp,Ki=ki,retained_SG_MW=_cost(ctx,s,e),
      converted_GFL_MW=sum(ctx.power)-_cost(ctx,s,e),evaluation=ev,
      history=DataFrame(history),accepted_steps=accepted,guard,trust_radius)
end

function _checkpoint(path,s,y,ctx,ev,d,phase)
    e,kp,ki=decode(ctx,s,y);cost=_cost(ctx,s,e)
    data=Dict{String,Any}("phase"=>phase,"support"=>s,"epsilon"=>e,
        "rho"=>[begin j=findfirst(==(i),s);j===nothing ? 1.0 : 1-e[j] end for i in 30:39],
        "Kp"=>kp,"Ki"=>ki,"retained_SG_MW"=>cost,"disturbance_MW"=>d,
        "alpha"=>ev.alpha,"beta_sampled"=>ev.beta,"Finf_Hz"=>ev.Finf,
        "Fpeak_Hz"=>ev.Fpeak,
        "model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a")
    open(path,"w") do io;TOML.print(io,data);end
end

function kkt_audit(ctx,cand;rocof_window=nothing,disturbance_MW=100.0)
    s=cand.support;y=cand.y;ev=cand.evaluation
    J=constraint_jacobian(ctx,s,y,ev;rocof_window,disturbance_MW,
        fd_step=2e-4,active_tol=2e-4).J
    active=findall(abs.(ev.g).<=2e-4)
    ns=length(s); c=zeros(length(y));for (j,b) in enumerate(s);c[j]=ctx.power[b-29]/1000;end
    if isempty(active)
        return (;primal=max(0.0,maximum(ev.g)),stationarity=norm(c,Inf),
            complementarity=0.0,dual_feasibility=true,LICQ=true,SOSC="NOT_TESTED",
            multipliers=Float64[],active_rows=active,jacobian_rank=0)
    end
    Ja=J[active,:]
    μ=pinv(transpose(Ja))*(-c)
    # Use nonnegative least squares by dropping negative constraint multipliers.
    free=collect(eachindex(active))
    while !isempty(free)
        muf=pinv(transpose(Ja[free,:]))*(-c)
        if minimum(muf)>=-1e-7;μ=zeros(length(active));μ[free].=max.(muf,0);break;end
        deleteat!(free,argmin(muf));
    end
    stat=norm(c+transpose(Ja)*μ,Inf)
    ranka=rank(Ja;atol=1e-9)
    (;primal=max(0.0,maximum(ev.g)),stationarity=stat,
      complementarity=maximum(abs.(μ.*ev.g[active])),dual_feasibility=all(μ.>=-1e-8),
      LICQ=ranka==length(active),SOSC="HESSIAN_NOT_CERTIFIED",
      multipliers=μ,active_rows=active,jacobian_rank=ranka)
end

end
