using CSV, DataFrames, TOML, SHA, LinearAlgebra
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(@__DIR__, "event_screens", isempty(ARGS) ? "batch" : splitext(basename(ARGS[1]))[1])
mkpath(OUT)
const P = include(joinpath(ROOT, "experiments", "nonlinear_codesign_20261001", "PDPhysicalReference.jl"))
include(joinpath(ROOT, "experiments", "bnd_expQ2B", "certified_search", "ConsistentEventObservation.jl"))
const SIGMA_REQ = 0.05
const F_LIMIT = 0.5
const R_LIMIT = 0.5
const V_MIN = 0.9
const V_MAX = 1.1
const ACTUATOR_SLACK_MIN = 0.002
BLAS.set_num_threads(1)

function event_network(nw, bus, delta_mw)
    mdl = P.PD39.PD39Model
    vertices, edges = mdl.copy_network_components(nw)
    defaults = get_defaults_dict(vertices[bus])
    pk = only(filter(s -> occursin("Pset", string(s)), collect(keys(defaults))))
    qk = only(filter(s -> occursin("Qset", string(s)), collect(keys(defaults))))
    p0, q0 = defaults[pk], defaults[qk]
    affect = (u, p, ctx) -> begin
        p[pk] = p0 - delta_mw / 100
        p[qk] = q0
    end
    set_callback!(vertices[bus], PresetTimeComponentCallback([1.0], ComponentAffect(affect, (), (pk, qk))))
    active = Network(vertices, edges)
    set_jac_prototype!(active)
    active
end

function main()
    length(ARGS) in (2,3) || error("usage: run_event_screen.jl candidate.toml event_id [uncertain_event_MW]")
    path = abspath(ARGS[1]); d = TOML.parsefile(path)
    rho = Float64.(d["rho"]); kp = Float64.(d["Kp"]); ki = Float64.(d["Ki"])
    length(rho) == length(kp) == length(ki) == 10 || error("expected ten rho/Kp/Ki entries")
    all(0 .< rho .< 1) || error("Q0 audit requires one fixed interior architecture")
    kp0, ki0 = 10pi, (10pi)^2 / 4
    all(0.25kp0 .<= kp .<= 4kp0) || error("Kp outside frozen gain bounds")
    all(0.25ki0 .<= ki .<= 4ki0) || error("Ki outside frozen gain bounds")

    sha = bytes2hex(sha256(read(path)))
    println("Q0_BUILD candidate_sha256=", sha); flush(stdout)
    base = P.frozen_baseline()
    nw = P.build_architecture(base, rho, kp, ki)
    state = P.trim_state(nw, base, rho, kp, ki)
    trim = P.residual_audit(nw, state)
    trim.maximum < 1e-8 || error("equilibrium residual failed: $(trim.maximum)")
    audit = P.direct_power_audit(state, base, rho)
    all(audit.bounds_status .== "PASS") || error("equilibrium component/guard audit failed")
    CSV.write(joinpath(OUT, "Q0_TRIM_POWER_AUDIT.csv"), audit)

    lin = linearize_network(state)
    lp = jacobian_eigenvals(lin)
    ig = argmin(abs.(lp)); abs(lp[ig]) < 1e-5 || error("rotational gauge mode not found")
    phys = lp[setdiff(eachindex(lp), [ig])]
    alpha = maximum(real.(phys))
    CSV.write(joinpath(OUT, "Q0_FULL_PHYSICAL_SPECTRUM.csv"), DataFrame(real=real.(phys), imag=imag.(phys)))

    cases = [("bus8_minus100", 8, -100.0), ("bus16_plus100", 16, 100.0),
             ("bus16_minus100", 16, -100.0), ("bus29_plus100", 29, 100.0),
             ("bus29_minus100", 29, -100.0)]
    cases = filter(c -> c[1] == ARGS[2], cases)
    length(cases) == 1 || error("unknown event id")
    if length(ARGS)==3
        c=only(cases)
        cases=[(c[1],c[2],parse(Float64,ARGS[3]))]
    end
    dt = 0.01; tspan = (0.0, 61.0); ts = collect(0.0:dt:61.0)
    rows = NamedTuple[]
    for (label, bus, delta) in cases
        println("Q0_EVENT_START ", label); flush(stdout)
        active = event_network(nw, bus, delta)
        started = time(); lastlog = Ref(started)
        watch = DiscreteCallback((u,t,int)->time()-lastlog[]>20, int->begin
            println("Q0_PROGRESS ", label, " t=", int.t, " wall_s=", round(time()-started,digits=1)); flush(stdout)
            lastlog[] = time()
            time()-started > 240 ? terminate!(int) : SciMLBase.u_modified!(int, false)
        end; save_positions=(false,false))
        sol = solve(ODEProblem(active, state, tspan), Rodas5P();
            callback=CallbackSet(get_callbacks(active), watch),
            initializealg=SciMLBase.NoInit(), saveat=dt,
            abstol=1e-9, reltol=1e-9, maxiters=1_000_000)
        complete = SciMLBase.successful_retcode(sol.retcode) && sol.t[end] >= 61-1e-8
        if !complete
            row = (;event=label,bus,delta_load_MW=delta,complete=false,retcode=string(sol.retcode),
                last_time_s=sol.t[end],F_peak_Hz=NaN,RoCoF_peak_Hz_s=NaN,Vmin_pu=NaN,Vmax_pu=NaN,
                max_GFL_current_ratio=NaN,min_GFL_Vdc_pu=NaN,max_GFL_Vdc_pu=NaN,
                min_SG_actuator_fraction_slack=NaN,pass=false,runtime_s=time()-started)
            push!(rows,row); CSV.write(joinpath(OUT,"Q0_EVENT_METRICS.csv"),DataFrame(rows))
            error("event $label incomplete; cannot accept Q0")
        end

        nb=39; phase=zeros(nb,length(ts)); volt=similar(phase); previous=zeros(nb); angle=zeros(nb)
        current_max=0.0; vdc_min=Inf; vdc_max=-Inf; actuator_slack=Inf
        lag=round(Int,0.5/dt); lag*dt == 0.5 || error("window does not align")
        for (j,t) in enumerate(ts)
            te = abs(t-1.0)<1e-12 ? t+1e-9 : t
            ss = abs(t-1.0)<1e-12 ? consistent_event_observation(sol,te) : NWState(sol,te)
            for b in 1:nb
                ur=Float64(ss[VIndex(b,:busbar₊u_r)]); ui=Float64(ss[VIndex(b,:busbar₊u_i)])
                raw=atan(ui,ur)
                if j==1; previous[b]=raw; angle[b]=raw
                else; angle[b]+=mod(raw-previous[b]+pi,2pi)-pi; previous[b]=raw; end
                phase[b,j]=angle[b]; volt[b,j]=hypot(ur,ui)
            end
            for b in 30:39
                i=b-29
                if rho[i]>0
                    ir=Float64(ss[VIndex(b,:gfl₊filter₊i_f_r)]); ii=Float64(ss[VIndex(b,:gfl₊filter₊i_f_i)])
                    i0=hypot(Float64(state[VIndex(b,:gfl₊filter₊i_f_r)]),Float64(state[VIndex(b,:gfl₊filter₊i_f_i)]))
                    current_max=max(current_max,hypot(ir,ii)/max(i0,eps()))
                    vdc=Float64(ss[VIndex(b,:gfl₊v_dc_state)]); vdc_min=min(vdc_min,vdc); vdc_max=max(vdc_max,vdc)
                end
                if rho[i]<1 && b!=39
                    xg=Float64(ss[VIndex(b,:ctrld_gen₊gov₊xg1)]); lo=Float64(ss[VIndex(b,:ctrld_gen₊gov₊V_min)]); hi=Float64(ss[VIndex(b,:ctrld_gen₊gov₊V_max)])
                    vr=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr)]); vlo=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr_min)]); vhi=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr_max)])
                    actuator_slack=min(actuator_slack,(xg-lo)/(hi-lo),(hi-xg)/(hi-lo),(vr-vlo)/(vhi-vlo),(vhi-vr)/(vhi-vlo))
                end
            end
        end
        freq=zeros(nb,length(ts)); roc=similar(freq)
        for j in eachindex(ts),b in 1:nb
            p1=j>lag ? phase[b,j-lag] : 0.0; p2=j>2lag ? phase[b,j-2lag] : 0.0
            freq[b,j]=(phase[b,j]-p1)/(2pi*0.5)
            roc[b,j]=(phase[b,j]-2p1+p2)/(2pi*0.5^2)
        end
        firstpost=findfirst(t->t>=1.0,ts)
        fpeak=maximum(abs,view(freq,:,firstpost:length(ts)))
        rpeak=maximum(abs,view(roc,:,firstpost:length(ts)))
        vmin=minimum(view(volt,:,firstpost:length(ts))); vmax=maximum(view(volt,:,firstpost:length(ts)))
        pass=alpha<=-SIGMA_REQ && fpeak<=F_LIMIT && rpeak<=R_LIMIT && vmin>=V_MIN && vmax<=V_MAX && actuator_slack>=ACTUATOR_SLACK_MIN
        row=(;event=label,bus,delta_load_MW=delta,complete=true,retcode=string(sol.retcode),last_time_s=sol.t[end],
            F_peak_Hz=fpeak,RoCoF_peak_Hz_s=rpeak,Vmin_pu=vmin,Vmax_pu=vmax,
            max_GFL_current_ratio=current_max,min_GFL_Vdc_pu=vdc_min,max_GFL_Vdc_pu=vdc_max,
            min_SG_actuator_fraction_slack=actuator_slack,pass,runtime_s=time()-started)
        push!(rows,row); CSV.write(joinpath(OUT,"Q0_EVENT_METRICS.csv"),DataFrame(rows))
        println("Q0_EVENT_DONE ",row); flush(stdout)
    end
    feasible=alpha<=-SIGMA_REQ && all(getproperty(r,:pass) for r in rows)
    weights=[base.rows[b].P for b in 30:39]
    result=Dict{String,Any}(
        "status"=>(feasible ? "SCREEN_PASS_ONLY" : "SCREEN_FAIL"),
        "candidate_sha256"=>sha,"rho"=>rho,"Kp"=>kp,"Ki"=>ki,
        "GFL_MW"=>dot(weights,rho),"GFL_percent"=>100*dot(weights,rho)/sum(weights),
        "retained_SG_MW"=>dot(weights,1 .-rho),"physical_eigenvalues"=>length(phys),
        "critical_real_part_s_inv"=>alpha,"equilibrium_residual_inf"=>trim.maximum,
        "all_five_events_pass"=>false,"screen_event_pass"=>all(getproperty(r,:pass) for r in rows),
        "gain_bounds_pass"=>true,"events"=>[Dict(string(k)=>getproperty(r,k) for k in propertynames(r)) for r in rows]
    )
    open(joinpath(OUT,"Q0_RESULT.toml"),"w") do io; TOML.print(io,result); end
    println("Q0_RESULT status=",result["status"]," retained_SG_MW=",result["retained_SG_MW"]," GFL_percent=",result["GFL_percent"])
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()


