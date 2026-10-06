module ReducedDAE

using LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff, SciMLBase, OrdinaryDiffEqRosenbrock
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
const N=PDExactDesignN
const SG=N.CollectiveModel.AnalyticSG
const OUT=joinpath(ROOT,"reports","nonlinear_codesign_20261001")

"Exact stator affine map; the machine's internal nonlinear dynamics remain intact."
function norton(x,p)
    pq,pd,ed,eq,w,delta=length(x)==12 ? view(x,7:12) : x
    ad=(p.xdpp-p.xls)/(p.xdp-p.xls); aq=(p.xqpp-p.xls)/(p.xqp-p.xls)
    cd=ad*eq+(1-ad)*pd; cq=-aq*ed+(1-aq)*pq
    sn,cs=sincos(delta); T=[sn -cs;cs sn]
    R=[p.rs -w*p.xqpp;w*p.xdpp p.rs]
    scale=p.rating_mva/p.system_base_mva
    -scale*T'*(R\T),scale*T'*(R\(w.*[-cq,cd]))
end

"Kron elimination of passive buses for a fixed Z-load disturbance. No dynamic truncation."
function network(ctx,bus=8,delta=0.0)
    Y=copy(ctx.net.y_static)
    loads=CSV.read(joinpath(ROOT,"reports","experiment_D","inputs","load.csv"),DataFrame)
    row=only(eachrow(loads[loads.bus.==bus,:]))
    audit=only(eachrow(ctx.net.load_audit[ctx.net.load_audit.bus.==bus,:]))
    bus in (31,39) && error("Generator-load buses require independently frozen initialized Pset")
    vset2=Float64(row.Pset)/Float64(audit.initialized_admittance_real)
    vset2>0 || error("invalid recovered Vset")
    for i in (2bus-1,2bus);Y[i,i]-=delta/(100vset2);end
    # All dynamic injections are at buses 30:39.
    passive=1:58; ports=59:78
    lift=-(Y[passive,passive]\Y[passive,ports])
    reduced=Y[ports,ports]+Y[ports,passive]*lift
    (;Y,reduced,lift,vset2,bus,delta)
end

function model(ctx,rho,kp,ki;bus=8,delta=0.,dc_convention=:legacy)
    all(0 .<=rho.<=1) || error("invalid replacement fraction")
    sgidx=UnitRange{Int}[]; gfidx=UnitRange{Int}[]; x0=Float64[]
    sp=[ctx.net.sg[b].op.parameters for b in 30:39]
    gp=[N.trim_gfl(ctx,b).pars for b in 30:39]
    for i in 1:10
        start=length(x0)+1
        rho[i]<1 && append!(x0,ctx.net.sg[i+29].op.x)
        push!(sgidx,start:length(x0))
        start=length(x0)+1
        rho[i]>0 && append!(x0,N.trim_gfl(ctx,i+29).x)
        push!(gfidx,start:length(x0))
    end
    dc_convention in (:legacy,:physical_supply) || error("unknown DC convention")
    (;ctx,rho=copy(rho),kp=copy(kp),ki=copy(ki),sgidx,gfidx,sp,gp,x0,net=network(ctx,bus,delta),dc_convention)
end

"Positive AC power flows from the converter to its filter; positive Pdc supplies it."
function gfl_rhs(x,u,p,Kp,Ki,convention)
    convention==:legacy && return N.gfl_rhs(x,u,p,Kp,Ki)
    gammaq,gammad,theta,omega,xi,ifi,ifr,vdi,vdc=x
    ur,ui=u;c=cos(theta);s=sin(theta)
    id=c*ifr+s*ifi;iq=-s*ifr+c*ifi
    error_dc=vdc-p.Vdc
    ed=error_dc*p.dc_kp+vdi-id;eq=p.iset_q-iq
    vid=p.cc_kp*ed+p.cc_ki*gammad;viq=p.cc_kp*eq+p.cc_ki*gammaq
    vir=c*vid-s*viq;vii=s*vid+c*viq;pac=vid*id+viq*iq;e=-s*ur+c*ui
    [eq,ed,omega,(xi+Kp*e-omega)/p.pll_tau,Ki*e,
        (p.omega_base/p.Xf)*(vii-ui-p.Rf*ifi-p.omega_frame*p.Xf*ifr),
        (p.omega_base/p.Xf)*(vir-ur-p.Rf*ifr+p.omega_frame*p.Xf*ifi),
        error_dc*p.dc_ki,(p.Pdc-pac)/(p.Cdc*vdc)]
end

function voltage(x,m;allbus=false)
    T=promote_type(eltype(x),eltype(m.rho),eltype(m.kp))
    G=Matrix{T}(m.net.reduced);h=zeros(T,20)
    for i in 1:10
        ix=(2i-1):(2i)
        if !isempty(m.sgidx[i])
            D,c=norton(view(x,m.sgidx[i]),m.sp[i])
            G[ix,ix].+=(1-m.rho[i]).*D;h[ix].+=(1-m.rho[i]).*c
        end
        if !isempty(m.gfidx[i])
            xi=m.gfidx[i];h[2i-1]+=m.rho[i]*x[xi[7]];h[2i]+=m.rho[i]*x[xi[6]]
        end
    end
    vp=-(G\h)
    allbus ? vcat(m.net.lift*vp,vp) : vp
end

function rhs!(dx,x,m,t)
    v=voltage(x,m)
    for i in 1:10
        u=view(v,(2i-1):(2i))
        if !isempty(m.sgidx[i])
            xi=m.sgidx[i];xs=view(x,xi);p=m.sp[i]
            dx[xi].=p.controlled ? SG.rhs(xs,u,p) : SG.rhs_uncontrolled(xs,u,p)
        end
        if !isempty(m.gfidx[i])
            xi=m.gfidx[i]
            dx[xi].=gfl_rhs(view(x,xi),u,m.gp[i],m.kp[i],m.ki[i],m.dc_convention)
        end
    end
    nothing
end

"Exact nonlinear Schur derivatives, in a fixed device support and limiter regime.
Parameters are (rho, log(Kp), log(Ki)); no time-domain finite differences are used.
"
function derivatives(x,m)
    nx=length(x);v=voltage(x,m)
    A=zeros(nx,nx);B=zeros(nx,20);C=zeros(20,nx);G=copy(m.net.reduced)
    gp=zeros(20,30);fp=zeros(nx,30)
    for i in 1:10
        vi=(2i-1):(2i);u=v[vi];sgcurrent=zeros(2);gfcurrent=zeros(2)
        if !isempty(m.sgidx[i])
            xi=m.sgidx[i];xx=x[xi];p=m.sp[i];n=length(xx)
            jf=ForwardDiff.jacobian(vcat(xx,u)) do z
                p.controlled ? SG.rhs(view(z,1:n),view(z,n+1:n+2),p) :
                    SG.rhs_uncontrolled(view(z,1:n),view(z,n+1:n+2),p)
            end
            jc=ForwardDiff.jacobian(z->SG.output(z,u,p),xx)
            A[xi,xi].=jf[:,1:n];B[xi,vi].=jf[:,n+1:n+2]
            C[vi,xi].+=(1-m.rho[i]).*jc
            D,c=norton(xx,p);G[vi,vi].+=(1-m.rho[i]).*D
            sgcurrent=SG.output(xx,u,p)
        end
        if !isempty(m.gfidx[i])
            xi=m.gfidx[i];xx=x[xi];p=m.gp[i]
            jf=ForwardDiff.jacobian(z->gfl_rhs(view(z,1:9),view(z,10:11),p,m.kp[i],m.ki[i],m.dc_convention),vcat(xx,u))
            A[xi,xi].=jf[:,1:9];B[xi,vi].=jf[:,10:11]
            C[2i-1,xi[7]]=m.rho[i];C[2i,xi[6]]=m.rho[i]
            gfcurrent=xx[[7,6]]
            e=-sin(xx[3])*u[1]+cos(xx[3])*u[2]
            fp[xi[4],10+i]=m.kp[i]*e/p.pll_tau
            fp[xi[5],20+i]=m.ki[i]*e
        end
        # Rho derivatives across support removal are not defined by this fixed-state map.
        if !isempty(m.sgidx[i]) && !isempty(m.gfidx[i])
            gp[vi,i]=gfcurrent-sgcurrent
        end
    end
    gv=G\hcat(C,gp)
    vx=-gv[:,1:nx];vp=-gv[:,nx+1:end]
    (;Fx=A+B*vx,Fp=fp+B*vp,v,vx,vp)
end

"Second-order L-stable SDIRK, with the exact derivative of its discrete trajectory.
Used to audit/compute search gradients. Final acceptance uses adaptive Rodas5P.
"
function tangent_simulate(m;horizon=60.,dt=.02,newton_tol=1e-10,sensitivities=true,monitor_buses=collect(30:39))
    all(0 .<m.rho.<1) || error("Tangent search requires an interior, fixed support")
    gamma=1-1/sqrt(2);h=dt;steps=round(Int,horizon/h);nx=length(m.x0)
    x=copy(m.x0);S=zeros(nx,30);xs=Vector{Vector{Float64}}(undef,steps+1)
    limiter_fraction=Inf;limiter_gradient=zeros(30)
    limiter_bus=0;limiter_kind="none";limiter_time=0.;limiter_side="none"
    vs=zeros(20,steps+1);dvs=Array{Float64}(undef,20,30,steps+1)
    xs[1]=copy(x);de=derivatives(x,m);vs[:,1]=de.v;dvs[:,:,1]=de.vp
    f=zeros(nx);a=h*gamma;b=h*(1-gamma)
    function stage(base,guess)
        y=copy(guess);laststep=Inf
        for it in 1:12
            rhs!(f,y,m,0.);rr=y-base-a*f
            de=derivatives(y,m);M=lu(I-a*de.Fx)
            norm(rr,Inf)<newton_tol && return y,copy(f),de,M
            dy=M\rr;y.-=dy;laststep=norm(dy,Inf)
        end
        error("SDIRK stage did not converge: last update $laststep")
    end
    for k in 1:steps
        y1,f1,d1,M1=stage(x,x)
        y2,f2,d2,M2=stage(x+b*f1,y1)
        if sensitivities
            T1=M1\(S+a*d1.Fp)
            S=M2\(S+b*(d1.Fx*T1+d1.Fp)+a*d2.Fp)
        end
        x=y2;xs[k+1]=copy(x);vs[:,k+1]=d2.v;dvs[:,:,k+1]=d2.vx*S+d2.vp
        for i in 1:10
            p=m.sp[i];p.controlled || continue;xi=m.sgidx[i]
            for (localidx,lo,hi) in ((1,p.gov_vmin,p.gov_vmax),(5,p.avr_vr_min,p.avr_vr_max))
                idx=xi[localidx];width=hi-lo
                for (slack,grad,side) in (((x[idx]-lo)/width,S[idx,:]/width,"lower"),((hi-x[idx])/width,-S[idx,:]/width,"upper"))
                    if slack<limiter_fraction
                        limiter_fraction=slack;limiter_gradient=copy(grad);limiter_bus=29+i
                        limiter_kind=localidx==1 ? "governor" : "AVR";limiter_time=k*h;limiter_side=side
                    end
                end
            end
        end
    end
    # Phase derivative is gauge-consistent; unwrapping adds parameter-independent integers.
    nb=length(monitor_buses);ph=zeros(nb,steps+1);dph=zeros(nb,30,steps+1)
    prev=angle.(m.ctx.net.voltage[monitor_buses]);init=copy(prev);unwrap=copy(prev)
    vmin=Inf;vmax=-Inf;gvmin=zeros(30);gvmax=zeros(30)
    lift=vcat(m.net.lift,Matrix{Float64}(I,20,20))
    for j in 1:steps+1
        vv=lift*vs[:,j];dv=lift*dvs[:,:,j]
        raw=atan.(vv[2 .*monitor_buses],vv[2 .*monitor_buses.-1]);unwrap.+=mod.(raw-prev.+pi,2pi).-pi
        ph[:,j]=unwrap-init;prev=raw
        for (i,bus) in enumerate(monitor_buses)
            ur=vv[2bus-1];ui=vv[2bus]
            dph[i,:,j]=(-ui*dv[2bus-1,:]+ur*dv[2bus,:])/(ur^2+ui^2)
        end
        for bus in 1:39
            ur=vv[2bus-1];ui=vv[2bus];vm=hypot(ur,ui)
            if vm<vmin;vmin=vm;gvmin=(ur*dv[2bus-1,:]+ui*dv[2bus,:])/vm;end
            if vm>vmax;vmax=vm;gvmax=(ur*dv[2bus-1,:]+ui*dv[2bus,:])/vm;end
        end
    end
    (;xs,vs,ph,dph,dt,times=collect(0:steps).*h,monitor_buses,vmin,vmax,gvmin,gvmax,
        limiter_fraction,limiter_gradient,limiter_bus,limiter_kind,limiter_time,limiter_side)
end

function tangent_peaks(r;window=.5)
    lag=round(Int,window/r.dt)
    abs(lag*r.dt-window)<1e-10 || error("window must be a multiple of the step")
    rows=NamedTuple[]
    for (name,coeff,denom) in (("F",(1.,-1.,0.),2pi*window),("R",(1.,-2.,1.),2pi*window^2))
        peak=-Inf;grad=zeros(30);bus=0;time=0.
        for j in axes(r.ph,2),i in axes(r.ph,1)
            val=0.;der=zeros(30)
            for q in 0:2
                jj=j-q*lag;jj>=1 || continue
                val+=coeff[q+1]*r.ph[i,jj]/denom
                der.+=coeff[q+1]*r.dph[i,:,jj]/denom
            end
            if abs(val)>peak
                peak=abs(val);grad=sign(val)*der;bus=r.monitor_buses[i];time=r.times[j]
            end
        end
        push!(rows,(;name,peak,gradient=grad,bus,time))
    end
    rows
end

function rotation_generator(x,m;jacobian=true)
    q=zeros(length(x));Dq=jacobian ? zeros(length(x),length(x)) : zeros(0,0)
    for i in 1:10
        !isempty(m.sgidx[i]) && (q[last(m.sgidx[i])]=1.)
        if !isempty(m.gfidx[i])
            ix=m.gfidx[i];q[ix[3]]=1.;q[ix[6]]=x[ix[7]];q[ix[7]]=-x[ix[6]]
            if jacobian;Dq[ix[6],ix[7]]=1.;Dq[ix[7],ix[6]]=-1.;end
        end
    end
    q,Dq
end

function rotate_state(x,m,angle)
    y=copy(x);s,c=sincos(angle)
    for i in 1:10
        !isempty(m.sgidx[i]) && (y[last(m.sgidx[i])]+=angle)
        if !isempty(m.gfidx[i])
            ix=m.gfidx[i];y[ix[3]]+=angle
            y[ix[7]]=c*x[ix[7]]-s*x[ix[6]];y[ix[6]]=s*x[ix[7]]+c*x[ix[6]]
        end
    end
    y
end

function simulate(m;horizon=60.,dt=.01,tol=1e-8,wall_limit=60.,maxiters=50000,rotating_frame=false)
    start=time()
    watch=DiscreteCallback((u,t,int)->time()-start>wall_limit,int->terminate!(int);save_positions=(false,false))
    # Use local device derivatives and the 20-port Schur solve instead of differentiating
    # a global 20x20 factorization in every state direction inside the stiff solver.
    ig=isempty(filter(!isempty,m.sgidx)) ? first(filter(!isempty,m.gfidx))[3] : last(last(filter(!isempty,m.sgidx)))
    xinit=copy(m.x0)
    if rotating_frame
        # The reference angle is replaced by the accumulated frame phase. This is
        # an invertible coordinate change, not a linearization or a truncated model.
        xinit[ig]=0.
        function canonical_rhs!(dz,z,m,t)
            y=copy(z);y[ig]=m.x0[ig];rhs!(dz,y,m,t)
            omega=dz[ig];q,_=rotation_generator(y,m;jacobian=false);dz.-=omega*q;dz[ig]=omega
        end
        function canonical_jac!(J,z,m,t)
            y=copy(z);y[ig]=m.x0[ig];f=similar(y);rhs!(f,y,m,t)
            omega=f[ig];A=derivatives(y,m).Fx;q,Dq=rotation_generator(y,m);ow=copy(A[ig,:])
            J.=A-q*ow'-omega*Dq;J[ig,:]=ow;J[:,ig].=0.
        end
        fun=ODEFunction(canonical_rhs!;jac=canonical_jac!)
    else
        fun=ODEFunction(rhs!;jac=(J,x,m,t)->(J.=derivatives(x,m).Fx))
    end
    elapsed=@elapsed raw=solve(ODEProblem(fun,xinit,(0.,horizon),m),Rodas5P();
        abstol=tol,reltol=tol,saveat=dt,dense=false,callback=watch,maxiters)
    ok=SciMLBase.successful_retcode(raw.retcode) && raw.t[end]>=horizon-1e-8
    if rotating_frame
        states=Vector{Vector{Float64}}(undef,length(raw.u))
        for j in eachindex(raw.u)
            y=copy(raw.u[j]);angle=y[ig];y[ig]=m.x0[ig];states[j]=rotate_state(y,m,angle)
        end
        sol=(;u=states,t=raw.t)
    else
        sol=raw
    end
    (;sol,ok,elapsed,retcode=string(raw.retcode),last_time=raw.t[end],rotating_frame)
end

function metrics(m,r;dt=.01,window=.5,savepath=nothing,monitor_buses=collect(30:39))
    n=length(r.sol.t);nb=length(monitor_buses);angles=zeros(nb,n);volt=zeros(39,n)
    currents=Float64[];dcs=Float64[];limiter_slack=Inf;limiter_fraction=Inf
    limiter_bus=0;limiter_kind="none";limiter_time=0.;limiter_side="none"
    previous=angle.(m.ctx.net.voltage[monitor_buses]);initial=copy(previous);unwrapped=copy(previous)
    for j in 1:n
        x=r.sol.u[j];v=voltage(x,m;allbus=true)
        volt[:,j]=hypot.(v[1:2:end],v[2:2:end])
        raw=atan.(v[2 .*monitor_buses],v[2 .*monitor_buses.-1])
        unwrapped.+=mod.(raw-previous.+pi,2pi).-pi
        angles[:,j]=unwrapped-initial;previous=raw
        for i in 1:10
            if !isempty(m.gfidx[i])
                xx=view(x,m.gfidx[i]);base=m.x0[m.gfidx[i]]
                push!(currents,hypot(xx[6],xx[7])/hypot(base[6],base[7]));push!(dcs,xx[9])
            end
            if !isempty(m.sgidx[i]) && m.sp[i].controlled
                xx=view(x,m.sgidx[i]);p=m.sp[i]
                limiter_slack=min(limiter_slack,xx[1]-p.gov_vmin,p.gov_vmax-xx[1],xx[5]-p.avr_vr_min,p.avr_vr_max-xx[5])
                for (idx,lo,hi,kind) in ((1,p.gov_vmin,p.gov_vmax,"governor"),(5,p.avr_vr_min,p.avr_vr_max,"AVR"))
                    for (val,side) in (((xx[idx]-lo)/(hi-lo),"lower"),((hi-xx[idx])/(hi-lo),"upper"))
                        if val<limiter_fraction
                            limiter_fraction=val;limiter_bus=29+i;limiter_kind=kind;limiter_time=r.sol.t[j];limiter_side=side
                        end
                    end
                end
            end
        end
    end
    lag=round(Int,window/dt);freq=zeros(nb,n);roc=zeros(nb,n)
    for j in 1:n,i in 1:nb
        a=angles[i,j];b=j>lag ? angles[i,j-lag] : 0.;c=j>2lag ? angles[i,j-2lag] : 0.
        freq[i,j]=(a-b)/(2pi*window);roc[i,j]=(a-2b+c)/(2pi*window^2)
    end
    fi=argmax(abs.(freq));ri=argmax(abs.(roc))
    result=(;complete=r.ok,retcode=r.retcode,last_time=r.last_time,runtime_s=r.elapsed,
        Fpeak_Hz=maximum(abs,freq),Rpeak_Hz_s=maximum(abs,roc),
        F_bus=monitor_buses[fi[1]],F_time=r.sol.t[fi[2]],R_bus=monitor_buses[ri[1]],R_time=r.sol.t[ri[2]],
        Vmin=minimum(volt),Vmax=maximum(volt),
        current_ratio_max=isempty(currents) ? 0. : maximum(currents),
        DC_min=isempty(dcs) ? 2.5 : minimum(dcs),DC_max=isempty(dcs) ? 2.5 : maximum(dcs),
        limiter_slack,limiter_fraction,limiter_bus,limiter_kind,limiter_time,limiter_side,
        retained_MW=dot(m.ctx.power,1 .-m.rho))
    if savepath!==nothing
        df=DataFrame(time_after_event_s=r.sol.t,Fmax_Hz=vec(maximum(abs.(freq),dims=1)),
            Rmax_Hz_s=vec(maximum(abs.(roc),dims=1)),Vmin_pu=vec(minimum(volt,dims=1)),
            Vmax_pu=vec(maximum(volt,dims=1)),Vbus8_pu=volt[8,:])
        for i in 1:nb;df[!,Symbol("f_bus$(monitor_buses[i])_Hz")]=freq[i,:];end
        CSV.write(savepath,df)
    end
    result
end
end
