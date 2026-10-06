using LinearAlgebra, CSV, DataFrames, TOML, SHA, Statistics
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE
const R=ReducedDAE
const ROOT=R.ROOT
const OUT=joinpath(@__DIR__,"baseline_reproduction")
mkpath(OUT);BLAS.set_num_threads(1)
say(x...)= (println(x...);flush(stdout))
ctx=R.N.design_context(ROOT)
function load_candidate(id,path)
    d=TOML.parsefile(path)
    (;id,rho=Float64.(d["rho"]),kp=Float64.(d["Kp"]),ki=Float64.(d["Ki"]),
      source=relpath(path,ROOT),source_sha256=bytes2hex(sha256(read(path))))
end
p0=(;id="baseline_reconstructed_from_protocol",rho=fill(.875,10),
    kp=fill(R.N.K0P,10),ki=fill(R.N.K0I,10),
    source="reports/analytic_iteration_20261001/refined/protocol.toml seed",source_sha256="")
joint=load_candidate("joint_final",joinpath(@__DIR__,"inputs","prior_joint_input.toml"))
fixed=load_candidate("fixed_gains_final",joinpath(@__DIR__,"inputs","prior_fixed_input.toml"))
candidates=[p0,joint,fixed]
rows=NamedTuple[]
poles=NamedTuple[]
for c in candidates
    spec=R.N.spectrum(ctx,c.rho,c.kp,c.ki)
    vals=eigvals(spec.quotient'*spec.model.Ared*spec.quotient)
    m=R.model(ctx,c.rho,c.kp,c.ki;bus=8,delta=0.0,dc_convention=:physical_supply); f=zeros(length(m.x0))
    R.rhs!(f,m.x0,m,0.0)
    retained=dot(ctx.power,1 .-c.rho);total=sum(ctx.power);gfl=total-retained
    push!(rows,(;candidate=c.id,source=c.source,source_sha256=c.source_sha256,
        rho_weighted_percent=100dot(ctx.power,c.rho)/total,retained_SG_MW=retained,
        GFL_MW=gfl,Kp_min=minimum(c.kp),Kp_max=maximum(c.kp),Ki_min=minimum(c.ki),Ki_max=maximum(c.ki),
        equilibrium_rhs_inf=norm(f,Inf),physical_spectral_abscissa=maximum(real.(vals)),
        finite_pole_count=length(vals),spectrum_method="ReducedDAE exact fixed-support ODE Jacobian + gauge quotient",
        reproduction_status="COMPUTED"))
    for (j,v) in enumerate(vals)
        push!(poles,(;candidate=c.id,pole_id=j,real=real(v),imag=imag(v)))
    end
end
CSV.write(joinpath(OUT,"TABLE_D00_REDUCED_CANDIDATES.csv"),DataFrame(rows))
CSV.write(joinpath(OUT,"TABLE_D00_REDUCED_POLES.csv"),DataFrame(poles))
event_set=[("design_bus8_plus100",8,100.0,p0),("joint_bus8_plus100",8,100.0,joint),
("joint_bus8_minus100",8,-100.0,joint),("joint_bus16_plus100",16,100.0,joint),
("joint_bus16_minus100",16,-100.0,joint),("joint_bus29_plus100",29,100.0,joint),
("joint_bus29_minus100",29,-100.0,joint),("fixed_bus8_plus100",8,100.0,fixed)]
event_rows=NamedTuple[]
for (label,bus,delta,c) in event_set
    say("D00_TDS_START ",label)
    m=R.model(ctx,c.rho,c.kp,c.ki;bus,delta,dc_convention=:physical_supply)
    run=R.simulate(m;horizon=60.0,dt=.01,tol=2e-9,wall_limit=900.0,maxiters=1000000,rotating_frame=true)
    if !run.ok
        push!(event_rows,(;label,candidate=c.id,bus,delta_load_MW=delta,complete=false,retcode=run.retcode,
          last_time_s=run.last_time,runtime_s=run.elapsed,Fpeak_Hz=missing,Rpeak_Hz_s=missing,Vmin_pu=missing,Vmax_pu=missing,
          current_ratio_max=missing,DC_min=missing,DC_max=missing,actuator_slack=missing))
        CSV.write(joinpath(OUT,"TABLE_D00_REDUCED_EVENTS.csv"),DataFrame(event_rows))
        error("D00 simulation incomplete for $label at $(run.last_time)")
    end
    met=R.metrics(m,run;dt=.01,window=.5,monitor_buses=collect(1:39))
    push!(event_rows,(;label,candidate=c.id,bus,delta_load_MW=delta,complete=true,retcode=run.retcode,
        last_time_s=run.last_time,runtime_s=run.elapsed,Fpeak_Hz=met.Fpeak_Hz,Rpeak_Hz_s=met.Rpeak_Hz_s,
        Vmin_pu=met.Vmin,Vmax_pu=met.Vmax,current_ratio_max=met.current_ratio_max,
        DC_min=met.DC_min,DC_max=met.DC_max,actuator_slack=met.limiter_fraction,
        Fbus=met.F_bus,Rbus=met.R_bus,limiter_bus=met.limiter_bus,limiter_kind=met.limiter_kind))
    CSV.write(joinpath(OUT,"TABLE_D00_REDUCED_EVENTS.csv"),DataFrame(event_rows))
    say("D00_TDS_DONE ",label," F=",met.Fpeak_Hz," R=",met.Rpeak_Hz_s," V=",met.Vmin,"..",met.Vmax)
end
open(joinpath(OUT,"D00_REDUCED_STATUS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"REDUCED_MODEL_REPRODUCTION_COMPLETE",
      "model"=>"ReducedDAE fixed-support nonlinear model; static algebraic Schur elimination; no controller delay",
      "all_TDS_complete"=>all(x.complete for x in event_rows),
      "events"=>length(event_rows),
      "warning"=>"This is not full PowerDynamics parity until compared with independent compiled physical model. No delay claims." ))
end



