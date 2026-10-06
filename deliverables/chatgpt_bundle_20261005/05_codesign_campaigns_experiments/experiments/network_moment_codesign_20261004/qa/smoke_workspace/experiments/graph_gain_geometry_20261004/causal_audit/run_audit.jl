"""Exact algebraic gate only; reads old trajectories and runs no new dynamics."""
const AUDIT_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(AUDIT_DIR, "..", "..", ".."))
include(joinpath(REPO_ROOT, "experiments", "nonlinear_codesign_20261001", "ReducedDAE.jl"))
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, SHA, Dates
const R = ReducedDAE
const EVENTS = [(8,-100.), (16,100.), (16,-100.), (29,100.), (29,-100.)]
const RHO_GRID = [0., .875, .95, .99, .999, 1.]
const WINDOW = .5
BLAS.set_num_threads(1)

function main()
    started = time()
    ctx = R.N.design_context(REPO_ROOT)
    kp = fill(.9*2pi*5,10)
    ki = fill((2pi*5)^2/4,10)
    altkp = fill(1.2*2pi*5,10)
    altki = fill(.7*(2pi*5)^2/4,10)
    vstar = zeros(78)
    vstar[1:2:end] = real.(ctx.net.voltage)
    vstar[2:2:end] = imag.(ctx.net.voltage)
    norton_D = [R.norton(ctx.net.sg[b].op.x,ctx.net.sg[b].op.parameters)[1] for b in 30:39]
    rows = NamedTuple[]
    for rho in RHO_GRID
        pre = R.model(ctx,fill(rho,10),kp,ki;bus=8,delta=0.,dc_convention=:physical_supply)
        trim_error = maximum(abs,R.voltage(pre.x0,pre;allbus=true)-vstar)
        for (bus,delta) in EVENTS
            m = R.model(ctx,fill(rho,10),kp,ki;bus,delta,dc_convention=:physical_supply)
            vv = R.voltage(m.x0,m;allbus=true)
            vc = complex.(vv[1:2:end],vv[2:2:end])
            phase = angle.(vc./ctx.net.voltage)
            F = maximum(abs,phase)/(2pi*WINDOW)
            Rmetric = F/WINDOW
            volt = abs.(vc)
            Ge = copy(m.net.Y)
            for i in 1:10
                ix = (2*(i+29)-1):(2*(i+29))
                Ge[ix,ix] .+= (1-rho).*norton_D[i]
            end
            deltaY = m.net.Y-ctx.net.y_static
            formula_voltage = vstar-Ge\(deltaY*vstar)
            identity_error = maximum(abs,formula_voltage-vv)
            alt = R.model(ctx,fill(rho,10),altkp,altki;bus,delta,dc_convention=:physical_supply)
            gain_error = maximum(abs,R.voltage(alt.x0,alt;allbus=true)-vv)
            @assert max(trim_error,identity_error,gain_error)<1e-10
            push!(rows,(common_rho=rho,bus=bus,delta_MW=delta,
                F_zero_Hz=F,R_zero_Hz_s=Rmetric,Vmin_zero_pu=minimum(volt),
                Vmax_zero_pu=maximum(volt),max_phase_rad=maximum(abs,phase),
                max_phase_bus=argmax(abs.(phase)),
                instantaneous_gate_pass=F<=.5 && Rmetric<=.5 && minimum(volt)>=.9 && maximum(volt)<=1.1,
                trim_voltage_error_inf=trim_error,jump_identity_error_inf=identity_error,
                gain_voltage_error_inf=gain_error,condition_G_event=cond(Ge)))
        end
    end
    df = DataFrame(rows)
    CSV.write(joinpath(AUDIT_DIR,"TABLE_01_INSTANTANEOUS_ALGEBRAIC.csv"),df)
    earlyrows = NamedTuple[]
    trajectory_files = String[]
    for (bus,delta) in EVENTS
        path = joinpath(REPO_ROOT,"experiments","graph_gsp_codesign_20261003","eval","baseline",
            "bus$(bus)_$(Int(delta))_trajectory.csv")
        push!(trajectory_files,path)
        tr = CSV.read(path,DataFrame)
        early = tr[(tr.time_after_event_s .>= 0) .& (tr.time_after_event_s .<= .04+1e-12),:]
        @assert nrow(early)==5
        push!(earlyrows,(bus=bus,delta_MW=delta,sample_count=nrow(early),
            last_time_s=maximum(early.time_after_event_s),
            F_sample_max_Hz=maximum(early.Fmax_Hz),R_sample_max_Hz_s=maximum(early.Rmax_Hz_s),
            V_sample_min_pu=minimum(early.Vmin_pu),V_sample_max_pu=maximum(early.Vmax_pu)))
    end
    CSV.write(joinpath(AUDIT_DIR,"TABLE_02_EXISTING_BASELINE_FIRST_40MS.csv"),DataFrame(earlyrows))
    files = [joinpath(AUDIT_DIR,"PROTOCOL.md"),@__FILE__,
        joinpath(REPO_ROOT,"experiments","nonlinear_codesign_20261001","ReducedDAE.jl"),
        joinpath(REPO_ROOT,"experiments","graph_gsp_codesign_20261003","DelayedEvents.jl"),
        joinpath(REPO_ROOT,"experiments","graph_gsp_codesign_20261003","baseline.toml"),
        joinpath(REPO_ROOT,"src","bnd_model_expN","PDExactDesignN.jl"),
        joinpath(REPO_ROOT,"src","bnd_design_e","CollectiveModel.jl"),
        joinpath(REPO_ROOT,"src","bnd_design","AnalyticSG.jl"),
        joinpath(REPO_ROOT,"src","bnd_design","AnalyticGFLPLL.jl"),
        joinpath(REPO_ROOT,"reports","experiment_N","TABLE_N01_original_operating_point.csv"),
        joinpath(REPO_ROOT,"reports","experiment_N","TABLE_N01_original_operating_point.csv.sha256"),
        joinpath(REPO_ROOT,"reports","experiment_A","matrices","bus33_baseline_equilibrium.csv")]
    append!(files,[joinpath(REPO_ROOT,"reports","experiment_D","inputs",f) for f in
        ["bus.csv","branch.csv","load.csv","machine.csv","avr.csv","gov.csv"]])
    append!(files,trajectory_files)
    for f in ["Project.toml","Manifest.toml"]
        path=joinpath(REPO_ROOT,f);isfile(path) && push!(files,path)
    end
    hashes = Dict(replace(relpath(f,REPO_ROOT),'\\'=>'/')=>bytes2hex(sha256(read(f))) for f in files)
    output_hashes = Dict(f=>bytes2hex(sha256(read(joinpath(AUDIT_DIR,f)))) for f in
        ["TABLE_01_INSTANTANEOUS_ALGEBRAIC.csv","TABLE_02_EXISTING_BASELINE_FIRST_40MS.csv"])
    endpoint = df[df.common_rho .== 1.,:]
    manifest = Dict("status"=>"NEGATIVE_INSTANTANEOUS_GATE_NUMERICAL",
        "created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),
        "runtime_seconds"=>time()-started,"common_rho_grid"=>RHO_GRID,
        "events"=>[Dict("bus"=>b,"delta_MW"=>d) for (b,d) in EVENTS],
        "window_seconds"=>WINDOW,"cases"=>nrow(df),"all_cases_pass"=>all(df.instantaneous_gate_pass),
        "endpoint_all_events_pass"=>all(endpoint.instantaneous_gate_pass),
        "max_jump_identity_error_inf"=>maximum(df.jump_identity_error_inf),
        "max_trim_voltage_error_inf"=>maximum(df.trim_voltage_error_inf),
        "max_gain_voltage_error_inf"=>maximum(df.gain_voltage_error_inf),
        "source_sha256"=>hashes,"output_sha256"=>output_hashes)
    open(joinpath(AUDIT_DIR,"MANIFEST.toml"),"w") do io
        TOML.print(io,manifest;sorted=true)
    end
    println("cases=",nrow(df)," pass=",sum(df.instantaneous_gate_pass)," no_new_trajectories=true")
    println("rho=1 worst F=",maximum(endpoint.F_zero_Hz)," R=",maximum(endpoint.R_zero_Hz_s),
        " voltage range=",minimum(endpoint.Vmin_zero_pu),",",maximum(endpoint.Vmax_zero_pu))
    println("parity max: trim=",maximum(df.trim_voltage_error_inf),
        " jump=",maximum(df.jump_identity_error_inf)," gain=",maximum(df.gain_voltage_error_inf))
end

main()
