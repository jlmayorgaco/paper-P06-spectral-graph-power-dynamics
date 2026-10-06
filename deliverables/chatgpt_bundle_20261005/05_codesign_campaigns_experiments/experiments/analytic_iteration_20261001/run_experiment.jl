const STARTED=time()
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, SHA, ForwardDiff
const R=ReducedDAE
const ROOT=R.ROOT
const OUT=joinpath(ROOT,"reports","analytic_iteration_20261001","refined")
mkpath(OUT);BLAS.set_num_threads(1)
const CTX=R.N.design_context(ROOT)
const COST=vcat(-CTX.power/sum(CTX.power),zeros(20))
const LO=vcat(fill(.01,10),log.(CTX.kpmin),log.(CTX.kimin))
const HI=vcat(fill(.99,10),log.(CTX.kpmax),log.(CTX.kimax))
const SCALE=vcat(ones(10),fill(6.,20))
const LOG=NamedTuple[]
say(args...)=(println(args...);flush(stdout))
write_toml(name,data)=open(io->TOML.print(io,data),joinpath(OUT,name),"w")
unpack(p)=(p[1:10],exp.(p[11:20]),exp.(p[21:30]))
make(p;bus=8,delta=100.)=R.model(CTX,unpack(p)...;bus,delta,dc_convention=:physical_supply)
capacity(p)=100dot(CTX.power,p[1:10])/sum(CTX.power)

# Primal active-set QP, adapted from the existing Q2B solver. Every resulting
# working-set solution is separately checked against the explicit Schur formula.
function qp(H,c,A,b,z)
    W=Int[];n=length(c)
    for it in 1:500
        if isempty(W)
            d=-(H\(H*z+c));mu=Float64[]
        else
            AW=A[W,:];Z=nullspace(AW;rtol=1e-12)
            d=-Z*((Z'*H*Z)\(Z'*(H*z+c)))
            mu=pinv(AW';rtol=1e-12)*(-(H*(z+d)+c))
        end
        if norm(d,Inf)<2e-9
            if isempty(W)||minimum(mu)>=-1e-8
                return (;z,mu,W,status="OK",iterations=it)
            end
            deleteat!(W,argmin(mu));continue
        end
        alpha=1.;block=0
        for i in eachindex(b)
            i in W && continue
            den=dot(A[i,:],d);den<=1e-10 && continue
            t=(b[i]-dot(A[i,:],z))/den
            if t<alpha;alpha=max(0.,t);block=i;end
        end
        z+=alpha*d
        if block!=0
            (isempty(W)||rank(A[vcat(W,block),:];rtol=1e-10)>length(W)) ||
                return (;z,mu,W,status="DEPENDENT",iterations=it)
            push!(W,block)
        end
    end
    (;z,mu=Float64[],W,status="ITERATION_LIMIT",iterations=500)
end

function modal(p;grad=false)
    rho,kp,ki=unpack(p);s=R.N.spectrum(CTX,rho,kp,ki)
    vals=eigen(s.quotient'*s.model.Ared*s.quotient).values
    ids=sort([i for i in eachindex(vals) if imag(vals[i])>=-1e-8];by=i->real(vals[i]),rev=true)[1:3]
    g=(real.(vals[ids]).+.05)/.05
    if grad
        rows=[begin
            d=R.N.simple_mode_sensitivities(CTX,rho,kp,ki;mode=i)
            d.condition<1e6 || error("Ill-conditioned mode")
            real.(vcat(d.rho,d.Kp.*kp,d.Ki.*ki))/.05
        end for i in ids]
        return (;g,J=reduce(vcat,transpose.(rows)),alpha=s.alpha)
    end
    (;g,alpha=s.alpha)
end

function evaluate(p,label;horizon=60.,tol=2e-9,save=true,bus=8,delta=100.,sample=.01)
    say("EVAL_START ",label," bus=",bus," delta=",delta)
    m=make(p;bus,delta)
    run=R.simulate(m;horizon,dt=sample,tol,wall_limit=150.,maxiters=500000,rotating_frame=true)
    run.ok || error("Incomplete trajectory $label: $(run.retcode), t=$(run.last_time)")
    met=R.metrics(m,run;dt=sample,monitor_buses=collect(1:39),
        savepath=save ? joinpath(OUT,label*"_trajectory.csv") : nothing)
    mo=modal(p)
    g=vcat(mo.g,[met.Fpeak_Hz/.5-1,met.Rpeak_Hz_s/.5-1,(.9-met.Vmin)/.1,
                 (met.Vmax-1.1)/.1,(.002-met.limiter_fraction)/.1])
    row=(;label,bus,delta,capacity_percent=capacity(p),alpha=mo.alpha,max_constraint=maximum(g),met...)
    push!(LOG,row);CSV.write(joinpath(OUT,"evaluations.csv"),DataFrame(LOG))
    say("EVAL_DONE ",label," C=",capacity(p)," F=",met.Fpeak_Hz," R=",met.Rpeak_Hz_s,
        " alpha=",mo.alpha," maxg=",maximum(g)," time=",run.elapsed)
    (;g,met,alpha=mo.alpha)
end

function gradient(p,label;horizon=20.,dt=.0125)
    say("TANGENT_START ",label," T=",horizon," dt=",dt)
    elapsed=@elapsed tr=R.tangent_simulate(make(p);horizon,dt,monitor_buses=collect(1:39))
    tr.limiter_fraction>0 || error("Smooth tangent invalid: actuator limiter reached")
    peaks=R.tangent_peaks(tr);mo=modal(p;grad=true)
    g=vcat(mo.g,[peaks[1].peak/.5-1,peaks[2].peak/.5-1,(.9-tr.vmin)/.1,
                 (tr.vmax-1.1)/.1,(.002-tr.limiter_fraction)/.1])
    J=vcat(mo.J,peaks[1].gradient'/.5,peaks[2].gradient'/.5,-tr.gvmin'/.1,tr.gvmax'/.1,-tr.limiter_gradient'/.1)
    CSV.write(joinpath(OUT,label*"_jacobian.csv"),DataFrame(J,:auto))
    write_toml(label*"_tangent.toml",Dict("constraints"=>g,"horizon_s"=>horizon,"step_s"=>dt,
        "runtime_s"=>elapsed,"F_time"=>peaks[1].time,"R_time"=>peaks[2].time,
        "F_bus"=>peaks[1].bus,"R_bus"=>peaks[2].bus))
    say("TANGENT_DONE ",label," F=",peaks[1].peak," R=",peaks[2].peak," time=",elapsed)
    (;g,J,tr,peaks,elapsed)
end

function gradient_gate(p)
    m=make(p);x=copy(m.x0);x.+=1e-6*sin.(eachindex(x));d=R.derivatives(x,m)
    auto=ForwardDiff.jacobian(p) do q
        mm=merge(m,(;rho=q[1:10],kp=exp.(q[11:20]),ki=exp.(q[21:30])))
        f=similar(q,length(x));R.rhs!(f,x,mm,0.);f
    end
    local_error=norm(d.Fp-auto,Inf)/max(norm(auto,Inf),eps())
    local_error<1e-10 || error("Local DAE derivative gate failed")
    # The initialized per-unit state/voltage is analytically independent of p
    # in this fixed-sharing trim; check its residual and DC sign similarity.
    nominal=make(p;delta=0.);f=similar(nominal.x0);R.rhs!(f,nominal.x0,nominal,0.)
    trim_error=norm(f,Inf);trim_error<1e-8 || error("Trim gate failed")
    Aphys=R.derivatives(nominal.x0,nominal).Fx
    legacy=R.N.descriptor(CTX,unpack(p)...).Ared
    signs=ones(length(nominal.x0));for ix in nominal.gfidx;signs[ix[9]]=-1.;end
    similarity=norm(Aphys-(signs.*legacy).*signs')/norm(Aphys)
    similarity<1e-10 || error("Physical modal similarity gate failed")
    tg=gradient(p,"gradient_gate";horizon=10.,dt=.0125)
    rows=NamedTuple[]
    for j in (10,20,30),h in (1e-4,5e-5)
        pp=copy(p);pm=copy(p);pp[j]+=h;pm[j]-=h
        ap=evaluate(pp,"fd_$(j)_$(h)_plus";horizon=10.,tol=2e-11,save=false,sample=.0125)
        am=evaluate(pm,"fd_$(j)_$(h)_minus";horizon=10.,tol=2e-11,save=false,sample=.0125)
        for (metric,k,vp,vm) in (("F",1,ap.met.Fpeak_Hz,am.met.Fpeak_Hz),("R",2,ap.met.Rpeak_Hz_s,am.met.Rpeak_Hz_s))
            fd=(vp-vm)/(2h);analytic=tg.peaks[k].gradient[j]
            rel=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-5)
            push!(rows,(;parameter=j,metric,h,analytic,finite_difference=fd,relative_error=rel))
        end
        CSV.write(joinpath(OUT,"gradient_gate.csv"),DataFrame(rows))
    end
    maxerr=maximum(r.relative_error for r in rows)
    write_toml("gradient_gate_summary.toml",Dict("local_parameter_derivative_error"=>local_error,
        "trim_residual"=>trim_error,"physical_modal_similarity_error"=>similarity,
        "trajectory_gradient_relative_error"=>maxerr,"threshold"=>.02,
        "scope"=>"Three parameter classes at bus39; nonlinear SDIRK tangents versus independently integrated adaptive Rodas5P finite differences; no global claim"))
    maxerr<.02 || error("Trajectory gradient gate failed: $maxerr; refine tangent step")
    say("GRADIENT_GATE_PASS error=",maxerr)
end

function save_candidate(p,label,ev)
    rho,kp,ki=unpack(p)
    write_toml(label*".toml",Dict("rho"=>rho,"Kp"=>kp,"Ki"=>ki,"dc_convention"=>"physical_supply",
        "replacement_percent"=>capacity(p),"retained_SG_MW"=>dot(CTX.power,1 .-rho),
        "constraints"=>ev.g,"status"=>"FINITE_EXPERIMENT_FEASIBLE_NOT_OPTIMUM"))
end

function optimize(p0,ev0,arm;iterations=4)
    p=copy(p0);ev=ev0;ids=arm=="joint" ? collect(1:30) : collect(1:10)
    scales=SCALE[ids];n=length(ids);trust=.015;history=NamedTuple[]
    push!(history,(;iteration=0,capacity_percent=capacity(p),retained_MW=ev.met.retained_MW,
        F=ev.met.Fpeak_Hz,R=ev.met.Rpeak_Hz_s,alpha=ev.alpha,maxg=maximum(ev.g),step=0.,explicit_error=0.,elapsed_s=0.))
    started=time()
    for k in 1:iterations
        tg=gradient(p,"$(arm)_$(k)")
        if any((ev.g.>-.1).&(abs.(ev.g-tg.g).>.004))
            tg=gradient(p,"$(arm)_$(k)_full";horizon=60.)
            any((ev.g.>-.1).&(abs.(ev.g-tg.g).>.004)) && error("Tangent refinement required")
        end
        A=tg.J[:,ids].*scales';H=Matrix{Float64}(I,n,n)*.01;c=COST[ids].*scales
        accepted=false;explicit_error=NaN
        for attempt in 1:4
            lower=max.((LO[ids]-p[ids])./scales,-trust);upper=min.((HI[ids]-p[ids])./scales,trust)
            Aq=vcat(A,Matrix{Float64}(I,n,n),-Matrix{Float64}(I,n,n))
            # Measured continuous-trajectory values anchor the affine constraints.
            bq=vcat(-ev.g.-.001,upper,-lower)
            minimum(bq)>=0 || error("Current iterate lacks the declared proposal guard")
            sol=qp(H,c,Aq,bq,zeros(n))
            if sol.status!="OK";say("QP_RETRY ",sol.status);trust*=.5;continue;end
            z=-(H\c)
            if !isempty(sol.W)
                AW=Aq[sol.W,:];M=AW*(H\AW')
                z.-=(H\AW')*(M\(-bq[sol.W]-AW*(H\c)))
            end
            explicit_error=norm(z-sol.z,Inf)
            explicit_error<1e-7 || error("Explicit iteration disagrees with QP: $explicit_error")
            maximum(Aq*sol.z-bq)<1e-7 || error("QP primal violation")
            step=zeros(30);step[ids]=scales.*sol.z
            -dot(COST,step)>1e-8 || break
            for bt in 0:5
                eta=2.0^(-bt);trial=p+eta*step
                tag="$(arm)_$(k)_$(attempt)_$(bt)"
                proposal=try
                    evaluate(trial,tag)
                catch err
                    startswith(sprint(showerror,err),"Incomplete trajectory") || rethrow()
                    write_toml(tag*"_incomplete.toml",Dict("error"=>sprint(showerror,err),"physical_instability_proven"=>false))
                    say("REJECT_INCOMPLETE ",tag);continue
                end
                # Guard is retained to make the next QP start strictly feasible.
                if maximum(proposal.g)<=-.00101
                    p=trial;ev=proposal;accepted=true
                    push!(history,(;iteration=k,capacity_percent=capacity(p),retained_MW=ev.met.retained_MW,
                        F=ev.met.Fpeak_Hz,R=ev.met.Rpeak_Hz_s,alpha=ev.alpha,maxg=maximum(ev.g),
                        step=norm(eta*step,Inf),explicit_error,elapsed_s=time()-started))
                    CSV.write(joinpath(OUT,arm*"_history.csv"),DataFrame(history))
                    save_candidate(p,arm*"_final",ev)
                    say("ACCEPT ",arm," k=",k," C=",capacity(p)," formula_error=",explicit_error)
                    break
                end
            end
            accepted && break
            trust*=.5
        end
        accepted || (say("STOP_NO_ACCEPTED_STEP ",arm);break)
    end
    CSV.write(joinpath(OUT,arm*"_history.csv"),DataFrame(history));save_candidate(p,arm*"_final",ev)
    p,ev
end

function main()
    paths=["Project.toml","Manifest.toml","experiments/analytic_iteration_20261001/run_experiment.jl",
        "experiments/nonlinear_codesign_20261001/ReducedDAE.jl","experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl",
        "src/bnd_model_expN/PDExactDesignN.jl","reports/experiment_Q2B/CERTIFIED_SEARCH/Z_LOCAL_SECURE_FINAL.toml"]
    hashes=Dict(path=>bytes2hex(sha256(read(joinpath(ROOT,path)))) for path in paths)
    write_toml("source_hashes_before.toml",hashes)
    write_toml("protocol.toml",Dict("network"=>"IEEE39","dc_convention"=>"physical_supply",
        "design_event_bus"=>8,"design_delta_MW"=>100.,"measurement_window_s"=>.5,
        "monitor_buses"=>collect(1:39),"validation_horizon_s"=>60.,"tangent_horizon_s"=>20.,"tangent_step_s"=>.0125,
        "frequency_limit_Hz"=>.5,"rocof_limit_Hz_s"=>.5,"voltage_min_pu"=>.9,"voltage_max_pu"=>1.1,
        "modal_decay_s_inv"=>.05,"minimum_actuator_fraction_slack"=>.002,
        "accepted_iterations_per_arm"=>4,"arms"=>["fixed_gains","joint"],
        "seed"=>"rho=0.875 at every generator; nominal PLL gains; original rho=0.9 failed preflight and is preserved in the parent directory",
        "continuous_disturbance_certificate"=>false,"current_and_energy_limits_certified"=>false,
        "global_or_local_optimality_claim"=>false))
    p0=vcat(fill(.875,10),fill(log(R.N.K0P),10),fill(log(R.N.K0I),10))
    say("BOOT_READY wall=",time()-STARTED)
    base=evaluate(p0,"baseline");maximum(base.g)<-.00101 || error("Declared baseline infeasible")
    save_candidate(p0,"baseline",base)
    gradient_gate(p0)
    joint,jv=optimize(p0,base,"joint")
    fixed,fv=optimize(p0,base,"fixed_gains")
    # Refined numerical acceptance and finite external-event falsification.
    evaluate(joint,"joint_refined";tol=2e-11,sample=.005)
    evaluate(fixed,"fixed_refined";tol=2e-11,sample=.005)
    for (b,d) in ((8,-100.),(16,100.),(16,-100.),(29,100.),(29,-100.))
        tag="holdout_$(b)_$(Int(d))"
        try
            evaluate(joint,tag;bus=b,delta=d)
        catch err
            startswith(sprint(showerror,err),"Incomplete trajectory") || rethrow()
            write_toml(tag*"_incomplete.toml",Dict("error"=>sprint(showerror,err),"physical_instability_proven"=>false))
            say("HOLDOUT_INCOMPLETE ",tag)
        end
    end
    after=Dict(path=>bytes2hex(sha256(read(joinpath(ROOT,path)))) for path in paths)
    hashes==after || error("A protected source/candidate changed during the experiment")
    write_toml("source_hashes_after.toml",after)
    write_toml("finished.toml",Dict("baseline_capacity_percent"=>capacity(p0),"joint_capacity_percent"=>capacity(joint),
        "fixed_capacity_percent"=>capacity(fixed),"joint_gain_MW"=>sum(CTX.power)*(capacity(joint)-capacity(p0))/100,
        "joint_minus_fixed_MW"=>sum(CTX.power)*(capacity(joint)-capacity(fixed))/100,
        "protected_hashes_match"=>true,"elapsed_s"=>time()-STARTED,
        "scope"=>"Finite-budget nonlinear iteration experiment; holdouts are falsification, not constraints or robustness certificates."))
    say("EXPERIMENT_FINISHED wall=",time()-STARTED)
end

try
    main()
catch err
    write_toml("failure.toml",Dict("error"=>sprint(showerror,err),"elapsed_s"=>time()-STARTED))
    rethrow()
end
