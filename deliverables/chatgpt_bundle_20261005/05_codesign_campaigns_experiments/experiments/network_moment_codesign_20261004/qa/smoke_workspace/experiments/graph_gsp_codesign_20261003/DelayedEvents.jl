module DelayedEvents
using LinearAlgebra, CSV, DataFrames, TOML, SciMLBase, OrdinaryDiffEqRosenbrock
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
const R=ReducedDAE
const CASES=[(8,-100.),(16,100.),(16,-100.),(29,100.),(29,-100.)]
BLAS.set_num_threads(1)

function detector(x,m,de=nothing)
    v=de===nothing ? R.voltage(x,m) : de.v
    e=[-sin(x[m.gfidx[i][3]])*v[2i-1]+cos(x[m.gfidx[i][3]])*v[2i] for i=1:10]
    de===nothing && return e
    C=zeros(10,length(x));P=zeros(10,30)
    for i=1:10
        ix=m.gfidx[i][3];s,c=sincos(x[ix])
        C[i,:].=-s.*de.vx[2i-1,:]+c.*de.vx[2i,:]
        C[i,ix]+=-c*v[2i-1]-s*v[2i]
        P[i,:].=-s.*de.vp[2i-1,:]+c.*de.vp[2i,:]
    end
    e,C,P
end

function injection(m)
    B=zeros(length(m.x0),10)
    for i=1:10
        B[m.gfidx[i][4],i]=m.kp[i]/m.gp[i].pll_tau;B[m.gfidx[i][5],i]=m.ki[i]
    end
    B
end

function tangent(m;tau=.04,h=.02,horizon=60.,sensitivities=true)
    lag=round(Int,tau/h);abs(lag*h-tau)<1e-12 || error("delay must align with step")
    n=length(m.x0);steps=round(Int,horizon/h);gamma=1-1/sqrt(2);a=h*gamma;b=h*(1-gamma)
    B=injection(m);x=copy(m.x0);S=zeros(n,30)
    eh1=zeros(10,steps);eh2=zeros(10,steps+1)
    ep1=zeros(10,30,steps);ep2=zeros(10,30,steps+1)
    de0=R.derivatives(x,m);eh2[:,1],_,ep2[:,:,1]=detector(x,m,de0)
    phase=zeros(39,steps+1);dphase=zeros(39,30,steps+1)
    prev=angle.(m.ctx.net.voltage);initial=copy(prev);uw=copy(prev)
    lift=vcat(m.net.lift,Matrix{Float64}(I,20,20))
    vmin=Inf;vmax=-Inf;gvmin=zeros(30);gvmax=zeros(30)
    slack=Inf;gslack=zeros(30);slackbus=0;slacktime=0.
    function observe(j,x,S,de)
        vv=lift*de.v;dv=lift*(de.vx*S+de.vp)
        raw=atan.(vv[2:2:end],vv[1:2:end]);uw.+=mod.(raw-prev.+pi,2pi).-pi
        phase[:,j].=uw-initial;prev.=raw
        for i=1:39
            ur,ui=vv[2i-1],vv[2i];vm=hypot(ur,ui)
            dphase[i,:,j].=(-ui.*dv[2i-1,:]+ur.*dv[2i,:])./(vm^2)
            grad=(ur.*dv[2i-1,:]+ui.*dv[2i,:])./vm
            if vm<vmin;vmin=vm;gvmin=copy(grad);end
            if vm>vmax;vmax=vm;gvmax=copy(grad);end
        end
        for i=1:10
            p=m.sp[i];p.controlled || continue;ix=m.sgidx[i]
            for (jj,lo,hi) in ((1,p.gov_vmin,p.gov_vmax),(5,p.avr_vr_min,p.avr_vr_max))
                stateidx=ix[jj];width=hi-lo
                for (val,grad) in (((x[stateidx]-lo)/width,S[stateidx,:]/width),((hi-x[stateidx])/width,-S[stateidx,:]/width))
                    if val<slack;slack=val;gslack=copy(grad);slackbus=i+29;slacktime=(j-1)*h;end
                end
            end
        end
    end
    observe(1,x,S,de0)
    function stage(base,guess,ed,edp)
        y=copy(guess);f=zeros(n)
        for it=1:12
            R.rhs!(f,y,m,0.);de=R.derivatives(y,m);e,C,P=detector(y,m,de)
            f.+=B*(ed-e);J=de.Fx-B*C;fac=lu(I-a*J)
            residual=y-base-a*f
            if norm(residual,Inf)<1e-10
                fp=de.Fp+B*(edp-P)
                for i=1:10
                    ix=m.gfidx[i];fp[ix[4],10+i]+=m.kp[i]/m.gp[i].pll_tau*(ed[i]-e[i])
                    fp[ix[5],20+i]+=m.ki[i]*(ed[i]-e[i])
                end
                return y,copy(f),J,fp,fac,de,e,C,P
            end
            y.-=fac\residual
        end
        error("delayed SDIRK stage failed")
    end
    zero_error=zeros(10);zero_derivative=zeros(10,30)
    for k=1:steps
        e1=k>lag ? eh1[:,k-lag] : zero_error; p1=k>lag ? ep1[:,:,k-lag] : zero_derivative
        # At t=tau the delayed detector has its first jump. Integrate the
        # preceding interval with the left limit; the next stage uses the
        # post-event history. Using the right limit at this endpoint creates
        # an O(h) early impulse in an otherwise second-order stage formula.
        e2=k>lag ? eh2[:,k-lag+1] : zero_error; p2=k>lag ? ep2[:,:,k-lag+1] : zero_derivative
        y1,f1,J1,fp1,F1,de1,e,C,P=stage(x,x,e1,p1)
        S1=sensitivities ? F1\(S+a*fp1) : S
        eh1[:,k].=e;ep1[:,:,k].=C*S1+P
        y2,f2,J2,fp2,F2,de2,e,C,P=stage(x+b*f1,y1,e2,p2)
        S=sensitivities ? F2\(S+b*(J1*S1+fp1)+a*fp2) : S
        x=y2;eh2[:,k+1].=e;ep2[:,:,k+1].=C*S+P
        observe(k+1,x,S,de2)
    end
    trj=(;ph=phase,dph=dphase,dt=h,times=collect(0:steps).*h,monitor_buses=collect(1:39))
    peaks=R.tangent_peaks(trj;window=.5)
    vals=[peaks[1].peak,peaks[2].peak,vmin,vmax,slack]
    G=reduce(vcat,transpose.([peaks[1].gradient,peaks[2].gradient,gvmin,gvmax,gslack]))
    (;vals,G,slackbus,slacktime)
end

function adaptive(m;tau=.04,horizon=60.,dtmax=.01,tol=1e-9)
    B=injection(m);n=length(m.x0);sols=Any[];ends=Float64[];xstar=copy(m.x0)
    function history(t)
        t<=0 && return xstar
        j=clamp(searchsortedfirst(ends,t-1e-12),1,length(sols))
        t<=last(ends)+1e-9 || error("future history")
        sols[j](clamp(t,sols[j].t[1],sols[j].t[end]))
    end
    function fun!(dx,x,p,t)
        R.rhs!(dx,x,m,t)
        if tau>0
            ed=t-tau < -1e-12 ? zeros(10) : detector(history(max(0,t-tau)),m)
            dx.+=B*(ed-detector(x,m))
        end
    end
    function jac!(J,x,p,t)
        de=R.derivatives(x,m);J.=de.Fx
        if tau>0;_,C,_=detector(x,m,de);J.-=B*C;end
    end
    fn=ODEFunction(fun!;jac=jac!);t0=0.;x=copy(xstar);started=time()
    while t0<horizon-1e-10
        t1=min(horizon,tau>0 ? t0+tau : horizon)
        sol=solve(ODEProblem(fn,x,(t0,t1)),Rodas5P();reltol=tol,abstol=tol,dtmax,dense=true,save_everystep=true,maxiters=100000)
        SciMLBase.successful_retcode(sol) || error("DDE integration $(sol.retcode)")
        push!(sols,sol);push!(ends,t1);x=copy(sol.u[end]);t0=t1
        time()-started<400 || error("event exceeded wall limit")
    end
    ts=collect(0:.01:horizon);states=[history(t) for t in ts]
    r=(;sol=(;t=ts,u=states),ok=true,elapsed=time()-started,retcode="METHOD_OF_STEPS_RODAS5P",last_time=horizon)
    met=R.metrics(m,r;dt=.01,window=.5,monitor_buses=collect(1:39))
    met,r
end

function run_cases(command,path;h=.02,ctx=R.N.design_context(R.ROOT))
    path=abspath(path);id=splitext(basename(path))[1]
    d=TOML.parsefile(path)
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    out=joinpath(@__DIR__,command,id);mkpath(out);rows=NamedTuple[];grads=Matrix{Float64}[]
    tau=Float64(get(d,"tau",.04))
    for (bus,delta) in CASES
        println("EVENT_START ",command," ",id," ",bus," ",delta);flush(stdout)
        m=R.model(ctx,rho,kp,ki;bus,delta,dc_convention=:physical_supply)
        if command=="grad"
            started=time();r=tangent(m;tau,h)
            push!(rows,(;bus,delta,F=r.vals[1],R=r.vals[2],Vmin=r.vals[3],Vmax=r.vals[4],slack=r.vals[5],runtime=time()-started))
            push!(grads,r.G);CSV.write(joinpath(out,"gradients.csv"),DataFrame(reduce(vcat,grads),:auto))
        else
            met,r=adaptive(m;tau)
            push!(rows,(;bus,delta,F=met.Fpeak_Hz,R=met.Rpeak_Hz_s,Vmin=met.Vmin,Vmax=met.Vmax,slack=met.limiter_fraction,
                runtime=met.runtime_s,pass=met.Fpeak_Hz<=.5 && met.Rpeak_Hz_s<=.5 && met.Vmin>=.9 && met.Vmax<=1.1 && met.limiter_fraction>=.002))
            R.metrics(m,r;dt=.01,window=.5,monitor_buses=collect(1:39),savepath=joinpath(out,"bus$(bus)_$(Int(delta))_trajectory.csv"))
        end
        CSV.write(joinpath(out,"events.csv"),DataFrame(rows))
        println("EVENT_DONE ",last(rows));flush(stdout)
    end
    rows
end

function main()
    run_cases(ARGS[1],ARGS[2];h=length(ARGS)>=3 ? parse(Float64,ARGS[3]) : .02)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__);DelayedEvents.main();end
