using LinearAlgebra, CSV, DataFrames, TOML
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE
const R=ReducedDAE; const OUT=joinpath(@__DIR__,"baseline_reproduction")
BLAS.set_num_threads(1)

function physical_zip_event_model(ctx,rho,kp,ki;bus=8,delta=100.0)
    m0=R.model(ctx,rho,kp,ki;bus=8,delta=0.0,dc_convention=:physical_supply)
    Y=copy(ctx.net.y_static)
    # OfficialIEEE39's ZIP CSV omits Vset; PowerDynamics' ZIPLoad default is 1 pu.
    # For its frozen pure-Z fractions, changing Pset by delta/100 changes G by that amount.
    for i in (2bus-1,2bus); Y[i,i]-=delta/100; end
    passive=1:58; ports=59:78
    lift=-(Y[passive,passive]\Y[passive,ports])
    reduced=Y[ports,ports]+Y[ports,passive]*lift
    net=(;Y,reduced,lift,vset2=1.0,bus,delta)
    merge(m0,(;net))
end

function readcandidate(path)
    d=TOML.parsefile(path)
    (rho=Float64.(d["rho"]),kp=Float64.(d["Kp"]),ki=Float64.(d["Ki"]))
end
ctx=R.N.design_context(R.ROOT)
baseline=(rho=fill(.875,10),kp=fill(R.N.K0P,10),ki=fill(R.N.K0I,10))
joint=readcandidate(joinpath(@__DIR__,"inputs","prior_joint_input.toml"))
rows=NamedTuple[]
for (name,c,pdlabel) in (("baseline",baseline,"baseline"),("joint_final",joint,"joint_final"))
    m=physical_zip_event_model(ctx,c.rho,c.kp,c.ki)
    run=R.simulate(m;horizon=60.0,dt=.005,tol=2e-9,wall_limit=900.0,maxiters=1000000,rotating_frame=true)
    run.ok || error("Reduced ZIP-contract event failed: $(run.retcode) at $(run.last_time)")
    met=R.metrics(m,run;dt=.005,window=.5,monitor_buses=collect(1:39))
    push!(rows,(;candidate=name,full_pd_label=pdlabel,complete=run.ok,F_reduced=met.Fpeak_Hz,
        R_reduced=met.Rpeak_Hz_s,Vmin_reduced=met.Vmin,Vmax_reduced=met.Vmax,
        Fbus_reduced=met.F_bus,Rbus_reduced=met.R_bus))
end
CSV.write(joinpath(OUT,"TABLE_D00_EVENT_CONTRACT_DIAGNOSIS.csv"),DataFrame(rows))
println(DataFrame(rows))
